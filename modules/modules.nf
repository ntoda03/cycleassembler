#!/usr/bin/env nextflow

nextflow.enable.dsl=2


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*               NGMALIGN align reads to genome with ngm                     */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process NGMALIGN {
    input:
        tuple val(pair_id), path(reads)
        path reference

    output:
        tuple val(pair_id), path('align.ngm.bam'),                  emit: ngmbam

    script:
    def read_in = params.single_end ? "-q $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    ngm -b -r $reference $read_in -o align.ngm.bam -t $task.cpus > ngm.log 2> ngm.err
    """
}


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*               EXTRACTBAM extract reads from a bam file                    */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////


process EXTRACTBAM {
    input:
        tuple val(pair_id), path(bam)

    output:
        tuple val(pair_id), path('*.fq.gz'),     emit: extractread

    script:
    if( params.single_end ){
        """
        source $projectDir/bin/functions.sh
        extract_bam_reads_se $bam $task.cpus
        """
    }
    else {
        """
        source $projectDir/bin/functions.sh
        extract_bam_reads $bam $task.cpus
        """
    }
}

///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*               COMPLEXITYFILTER remove low complexity reads PE             */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process COMPLEXITYFILTER {
    input:
        tuple val(pair_id), path(reads)
        val trimargs

    output:
        tuple val(pair_id), path('good*fq.gz'),    emit: filterread

    script:
    def read_in = params.single_end ? "-i $reads" : "-i ${reads[0]} -I ${reads[1]}"
    def read_out = params.single_end ? "-o good.fq.gz" : "-o good.1.fq.gz -O good.2.fq.gz"
    """
    fastp -G -A -L -Q --low_complexity_filter $read_in $read_out
    """
}



///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*                       SPADESASSEM spades assembly                         */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process SPADESASSEM {
    input:
        tuple val(pair_id), path(reads)
        val assemargs

    output:
        tuple val(pair_id), path('spades.fa'),           emit: assembly

    script:
    def read_in = params.single_end ? "-s $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    spades.py $assemargs $read_in -t $task.cpus -o spades_assembly 
    if [ -s spades_assembly/scaffolds.fasta ]; then
      cp spades_assembly/scaffolds.fasta spades.fa
      echo "Assembly complete."
    else
      if [ -s spades_assembly/contigs.fasta ]; then
        cp spades_assembly/contigs.fasta spades.fa
        echo "Warning: Partial assembly, likely low quality data. Assembly complete."
      else
        echo "Error with assembly $pair_id."
        return 1
      fi
    fi
    """
}

///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*             BLASTFILTER only keep contigs that blast to reference         */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process BLASTFILTER {
    input:
        tuple val(pair_id), path(contigs)
        path reference
        val fasta_command

    output:
        tuple val(pair_id), path('scaffolds.verified.fasta'),           emit: filtercontigs

    script:
    """
    source $projectDir/bin/functions.sh
    $fasta_command -E 1e-10 -T $task.cpus -m 8  $contigs $reference | \
        sed 's/_/ /g' |sed 's/ /_/g' |awk '{print \$1}' |sort |uniq > scaffolds.blast.list
    if [ ! -s scaffolds.blast.list ]; then
      return 1
    fi
    extractBlastedScaff scaffolds.blast.list $contigs scaffolds.verified.fasta T
    """
}


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*             CYCLEASSEM iterative assembly based on updated reference      */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process CYCLEASSEM {
    publishDir "$params.outdir/assembled_contigs", mode: 'copy'

    input:
        tuple val(pair_id), path(initial_contigs)
        tuple val(pair_id), path(reads)
        path reference
        val fasta_command
        val maxit

    output:
        tuple val(pair_id), path("${pair_id}.final_scaffolds.fa"),           emit: cyclecontigs

    script:
    def ngm_in = params.single_end ? "-q $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    filter_in = params.single_end ? "-i reads/output.fq.gz" : "-i reads/output.1.fq.gz -I reads/output.2.fq.gz"
    filter_out = params.single_end ? "-o reads/good.fq.gz" : "-o reads/good.1.fq.gz -O reads/good.2.fq.gz"
    def spades_in = params.single_end ? "-s reads/good.fq.gz" : "-1 reads/good.1.fq.gz -2 reads/good.2.fq.gz"
    def bam_extract = params.single_end ? "extract_bam_reads_se" : "extract_bam_reads"
    """
    source $projectDir/bin/functions.sh
    i=0
    contcount=10 
    cycle_genome=$initial_contigs

    while ( [ \"\$i\" -lt \"$maxit\" ] && [ \"\$contcount\" -gt \"1\" ] ); do
        let i=i+1 
        mkdir run_\$i/

        #### Map reads against seed contigs to get mapped reads ####
        ngm -b -i 0.99 -r \$cycle_genome $ngm_in -o run_\$i/output.bam -t $task.cpus > run_\$i/ngm.log 2> run_\$i/ngm.err
        command_success=0
        grep '(0 reads mapped' run_\$i/ngm.err > /dev/null 2>1 || command_success=1
        if [ \"\$command_success\" -eq 0 ]; then
          if [ \"\$i\" -eq 1 ]; then
            echo 'Mapping to cycle assembly failed. Using inital contigs.'
            cp $initial_contigs ${pair_id}.final_scaffolds.fa
          else
            echo 'Subsequent cycle mapping failed. Using previous cycle.'
            let i=i-1
            cp run_\$i/scaffolds.verified.fasta ${pair_id}.final_scaffolds.fa
            break
          fi
        fi
        $bam_extract run_\$i/output $task.cpus
        mv run_\$i/*.fq.gz reads/
        rm -f run_\$i/*.bam run_\$i/*.ngm run_\$i/*.bt2

        #### Do de novo assembly of plastid reads ####
        fastp -G -A -L -Q --low_complexity_filter $filter_in $filter_out
        spades.py --cov-cutoff 1 $spades_in -t $task.cpus -o run_\$i/spades_assembly
        if [ -s run_\$i/spades_assembly/scaffolds.fasta ]; then
          cp run_\$i/spades_assembly/scaffolds.fasta run_\$i/scaff.fa
        else
          if [ -s run_\$i/spades_assembly/contigs.fasta ]; then
            echo 'Warning: only contigs produced in spades assembly.'
            cp run_\$i/spades_assembly/contigs.fasta run_\$i/scaff.fa
          else
            echo 'Error in spades cycle assembly. Low quality data likely.'
            exit 1
          fi
        fi

        #### Only keep scaffolds that blast back to reference genome to discard junk ####
        $fasta_command -E 1e-10 -T $task.cpus -m 8 run_\$i/scaff.fa $reference | \
            sed 's/_/ /g' |sed 's/ /_/g' |awk '{print \$1}' |sort |uniq > run_\$i/spades_assembly/scaffolds.blast.list
        if [ ! -s run_\$i/spades_assembly/scaffolds.blast.list ]; then
          if [ \"\$i\" -eq 1 ]; then
            echo 'No hits found. Using inital contigs.'
            cp $initial_contigs run_\$i/scaff.fa
          else
            echo 'No hits found. Using previous cycle.'
            let i=i-1
            cp run_\$i/scaffolds.verified.fasta run_\$i/scaff.fa
          fi
          $fasta_command -E 1e-10 -T $task.cpus -m 8 run_\$i/scaff.fa $reference | \
            sed 's/_/ /g' |sed 's/ /_/g' |awk '{print \$1}' |sort |uniq > run_\$i/spades_assembly/scaffolds.blast.list
          extractBlastedScaff run_\$i/spades_assembly/scaffolds.blast.list run_\$i/scaff.fa run_\$i/scaffolds.verified.fasta T
          samtools faidx run_\$i/scaffolds.verified.fasta
          break
        fi
        extractBlastedScaff run_\$i/spades_assembly/scaffolds.blast.list run_\$i/scaff.fa run_\$i/scaffolds.verified.fasta T

        #### Check summary statistics of new plastid assembly ####
        contcount=\$(grep -c \">\" run_\$i/scaffolds.verified.fasta)
        samtools faidx run_\$i/scaffolds.verified.fasta
        cycle_genome=run_\$i/scaffolds.verified.fasta
        rm -f run_\$i/*fq* run_\$i/*fa run_\$i/*ngm reads/*fq.gz
    done
    if [ \"\$i\" -eq \"$maxit\" ] || [ \"\$contcount\" -eq \"1\" ]; then
        cp run_\$i/scaffolds.verified.fasta ${pair_id}.final_scaffolds.fa
    fi
    """
}


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*            EXTRACTEXONS extract exon sequences based on reference         */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process EXTRACTEXONS {
    publishDir "$params.outdir/exons/sequences/", mode: 'copy'

    input:
        tuple val(pair_id), path(contigs)
        path(exons)

    output:
        tuple val(pair_id), path('${pair_id}.fa'),           emit: exonseqs

    script:
    """
    source $projectDir/bin/functions.sh
    cp $exons exons.fa
    # Extract sequences that are hits for exons
    sed -i 's/_length.*//g' $contigs
    # Coordinates start from the end of the contig if they align on negative strand so avoid that by adding revcomp of all contigs
    seqtk seq -r $contigs > scaffolds.rev.fa
    sed -i 's/>/>R/g' scaffolds.rev.fa
    cat scaffolds.rev.fa >> $contigs
    fasta36 -E 1E-10 -T $task.cpus -m 8 $exons $contigs > gene_search.txt 2> /dev/null
    awk '(\$8>\$7) && (\$10>\$9) {print \$0}' gene_search.txt > gene_search.stranded.txt
    # Only except hits that span >80% of the exon
    samtools faidx exons.fa
    $projectDir/bin/mergeBlastHits.py gene_search.stranded.txt gene_search.stranded_merge.txt flanking_positive 50
    join <(awk '{printf \"%s %s:%s-%s %s\\n\",\$1,\$2,\$9,\$10,\$8-\$7}' gene_search.stranded_merge.txt |sort -k1,1 |uniq) <(awk '{print \$1,\$2}' exons.fa.fai |sort -k1,1) \
        |awk '((0.8*\$4) < \$3) {print \$2}' |uniq > gene_search.filtered.txt
    extract_seq gene_search.filtered.txt $contigs ${pair_id}.fa F    
    sed -i 's/-/__/g' ${pair_id}.fa
    sed -i 's/:/___/g' ${pair_id}.fa
    """
}


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*             FINDEXONS find exon sequences in contigs                      */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process FINDEXONS {
    publishDir "$params.outdir/exons/identified/", mode: 'copy'

    input:
        tuple val(pair_id), path(seqs)
        path(exons)

    output:
        tuple val(pair_id), path("$pair_id/*.fa"),           emit: exonsbra

    script:
    """
    source $projectDir/bin/functions.sh
    # Get the best reciprical alignment between exons and extracted sequences to only have 1 per exons
    getBRA $seqs $exons dna_dna
    mkdir -p $pair_id/
    while IFS=' ' read col1 col2
    do
        echo \$col2
      samtools faidx $seqs \$col1 > $pair_id/\${col2}.fa
      sed -i 's/___/:/g' $pair_id/\${col2}.fa
      sed -i 's/__/-/g' $pair_id/\${col2}.fa
    done < ${seqs}.BRA.dna_dna.txt
    """
}


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*                 CLUSTER create clusters of exons                          */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process CLUSTER {
    publishDir "$params.outdir/exons/", mode: 'copy'

    input:
        tuple val(pair_id), path(seqs)
        path exons

    output:
        tuple val(pair_id), path("exon_clusters/$pair_id/*txt.fa"),           emit: clusterseqs
        tuple val(pair_id), path("exon_clusters/$pair_id/*clw"),              emit: alignments

    script:
    """
    source $projectDir/bin/functions.sh
    mkdir -p exon_clusters
    create_cluster $exons $seqs nucl exon_clusters/$pair_id $task.cpus
    cp exon_clusters/$pair_id/gene_search_clusters/*clw exon_clusters/$pair_id/
    cp exon_clusters/$pair_id/gene_search_clusters/*txt.fa exon_clusters/$pair_id/
   """
}