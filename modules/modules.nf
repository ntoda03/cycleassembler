#!/usr/bin/env nextflow

nextflow.enable.dsl=2


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*               NGMALIGN align reads to genome with ngm                     */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process NGMALIGN {
    input:
        tuple val(pair_id), path(reads1), path(reads2)
        path reference

    output:
        tuple val(pair_id), path('align.ngm.bam'),                  emit: ngmbam

    script:
    """
    ngm -b -r $reference -1 $reads1 -2 $reads2 -o align.ngm.bam -t $task.cpus > ngm.log 2> ngm.err
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
        tuple val(pair_id), path('align.ngm.1.fq'), path('align.ngm.2.fq'),     emit: extractread

    script:
    """
    source $projectDir/bin/functions.sh
    extract_bam_reads $bam $task.cpus
    """
}

///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*               COMPLEXITYFILTER remove low complexity reads PE             */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process COMPLEXITYFILTER {
    input:
        tuple val(pair_id), path(reads1), path(reads2)
        val trimargs

    output:
        tuple val(pair_id), path('good.1.fq.gz'), path('good.2.fq.gz'),     emit: filterread

    script:
    """
    fastp -G -A -L -Q --low_complexity_filter -i $reads1 -I $reads2 -o good.1.fq.gz -O good.2.fq.gz
    """
}



///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*                       SPADESASSEM spades assembly                         */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process SPADESASSEM {
    input:
        tuple val(pair_id), path(reads1), path(reads2)
        val assemargs

    output:
        tuple val(pair_id), path('spades.fa'),           emit: assembly

    script:
    """
    spades.py $assemargs -1 $reads1 -2 $reads2 -t $task.cpus -o spades_assembly 
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
    publishDir '$outdir/assmbled_contigs', mode: 'copy'

    input:
        tuple val(pair_id), path(initial_contigs)
        tuple val(pair_id), path(reads1), path(reads2)
        path reference
        val fasta_command
        val maxit

    output:
        tuple val(pair_id), path('final_scaffolds.fa'),           emit: cyclecontigs

    script:
    """
    source $projectDir/bin/functions.sh
    i=0
    contcount=10 
    cycle_genome=$initial_contigs

    while ( [ \"\$i\" -lt \"$maxit\" ] && [ \"\$contcount\" -gt \"1\" ] ); do
        let i=i+1 
        mkdir run_\$i/

        #### Map reads against seed contigs to get mapped reads ####
        ngm -b -i 0.99 -r \$cycle_genome -1 $reads1 -2 $reads2 -o run_\$i/output.bam -t $task.cpus > run_\$i/ngm.log 2> run_\$i/ngm.err
        command_success=0
        grep '(0 reads mapped' run_\$i/ngm.err > /dev/null 2>1 || command_success=1
        if [ \"\$command_success\" -eq 0 ]; then
          if [ \"\$i\" -eq 1 ]; then
            echo 'Mapping to cycle assembly failed. Using inital contigs.'
            cp $initial_contigs final_scaffolds.fa
          else
            echo 'Subsequent cycle mapping failed. Using previous cycle.'
            let i=i-1
            cp run_\$i/scaffolds.verified.fasta final_scaffolds.fa
            break
          fi
        fi
        extract_bam_reads run_\$i/output $task.cpus
        gzip run_\$i/output.*.fq
        rm -f run_\$i/*.bam run_\$i/*.ngm run_\$i/*.bt2

        #### Do de novo assembly of plastid reads ####
        fastp -G -A -L -Q --low_complexity_filter -i run_\$i/output.1.fq.gz -I run_\$i/output.2.fq.gz -o run_\$i/good.1.fq.gz -O run_\$i/good.2.fq.gz
        spades.py --cov-cutoff 1 -1 run_\$i/good.1.fq.gz -2 run_\$i/good.2.fq.gz -t $task.cpus -o run_\$i/spades_assembly
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
        rm -f run_\$i/*fq* run_\$i/*fa run_\$i/*ngm
    done
    if [ \"\$i\" -eq \"$maxit\" ] || [ \"\$contcount\" -eq \"1\" ]; then
        cp run_\$i/scaffolds.verified.fasta final_scaffolds.fa
    fi
    """
}


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*             MUSCLE align sequences with muscle      */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process MUSCLE {
    //input:
        //tuple val(pair_id), path(initial_contigs)

    //output:
        //tuple val(pair_id), path('final_scaffolds.fa'),           emit: cyclecontigs

    script:
    """
    muscle
    """
}


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*            EXTRACTEXONS extract exon sequences based on reference         */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process EXTRACTEXONS {
    input:
        tuple val(pair_id), path(contigs)
        path(exons)

    output:
        tuple val(pair_id), path('gene_search.fa'),           emit: exonseqs

    script:
    """
    cp $exon exons.fa
    # Extract sequences that are hits for exons
    sed -i 's/_length.*//g' $contigs
    # Coordinates start from the end of the contig if they align on negative strand so avoid that by adding revcomp of all contigs
    seqtk seq -r $contigs > scaffolds.rev.fa
    sed -i 's/>/>R/g' scaffolds.rev.fa
    cat scaffolds.rev.fa >> $contigs
    fasta36 -E 1E-10 -T $taks.cpus -m 8 $exons $contigs > gene_search.txt 2> /dev/null
    awk '(\$8>\$7) && (\$10>\$9) {print \$0}' gene_search.txt > gene_search.stranded.txt
    # Only except hits that span >80% of the exon
    samtools faidx exons.fa
    $projectDir/bin/mergeBlastHits.py gene_search.stranded.txt gene_search.stranded_merge.txt flanking_positive 50
    join <(awk '{printf \"%s %s:%s-%s %s\n\",\$1,\$2,\$9,\$10,\$8-\$7}' gene_search.stranded_merge.txt |sort -k1,1 |uniq) <(awk '{print \$1,\$2}' exons.fa.fai |sort -k1,1) \
        |awk '((0.8*\$4) < \$3) {print \$2}' |uniq > gene_search.filtered.txt
    extract_seq gene_search.filtered.txt $contigs gene_search.fa F    
    sed -i 's/-/__/g' gene_search.fa
    sed -i 's/:/___/g' gene_search.fa
    """
}


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*             FINDEXONS find exon sequences in contigs                      */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process FINDEXONS {
    publishDir '$outdir/exons', mode: 'copy', pattern: 'exon_sequences/*.fa'

    input:
        tuple val(pair_id), path(seqs)
        path(exons)

    output:
        tuple val(pair_id), path('exon_sequences/*.fa'),           emit: exonsbra

    script:
    """
    source $projectDir/bin/functions.sh
    # Get the best reciprical alignment between exons and extracted sequences to only have 1 per exons
    getBRA $seqs $exons dna_dna
    mkdir -p exon_sequences 
    while IFS=' ' read col1 col2
    do
      samtools faidx $seqs \$col1 > exon_sequences/\${col2}.fa
      sed -i 's/___/:/g' exon_sequences/\${col2}.fa
      sed -i 's/__/-/g' exon_sequences/\${col2}.fa
    done < ${seqs}.BRA.dna_dna.txt
    """
}

