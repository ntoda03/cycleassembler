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
    ngm -b -r $reference -1 $reads1 -2 $reads2 -o align.ngm.bam -t $task.cpus
    """
}


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*               EXTRACTBAM extract reads from a bam file                    */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////


process EXTRACTBAM {
    //beforeScript "source $projectDir/bin/functions.sh"

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

    //beforeScript "source $projectDir/bin/functions.sh"

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


