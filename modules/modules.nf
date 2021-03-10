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
        tuple val(pair_id), path('align.ngm.1.fq'), path('align.ngm.2.fq'),     emit: ngmread
        path 'align.ngm.bam',                                                     emit: ngmbam

    //beforeScript "source $projectDir/bin/functions.sh"

    script:
    """
    source $projectDir/bin/functions.sh
    ngm -b -r $reference -1 $reads1 -2 $reads2 -o align.ngm.bam -t $task.cpus
    extract_bam_reads align.ngm $task.cpus
    if [ -s align.ngm.1.fq ]; then
      echo 'Mapping complete.''
    else
      echo 'Error with mapping.''
      return 1  
    fi
    """
}


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*               AFTERQC automatic read filtering and QC                     */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process AFTERQC {
    input:
        tuple val(pair_id), path(reads1), path(reads2)
        val trimargs

    output:
        tuple val(pair_id), path('goodreads/*.1.good.fq'), path('goodreads/*.2.good.fq'),     emit: qcread

    script:
    """
    after.py -1 $reads1 -2 $reads2 $trimargs -g goodreads/ -b badreads/
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

    beforeScript "source $projectDir/bin/functions.sh"

    script:
    """
    $fasta_command -E 1e-10 -T $task.cpus -m 8  $contigs $reference | \
        sed 's/_/ /g' |sed 's/ /_/g' |awk '{print \$1}' |sort |uniq > scaffolds.blast.list
    if [ ! -s scaffolds.blast.list ]; then
      return 1
    fi
    extractBlastedScaff scaffolds.blast.list $contigs scaffolds.verified.fasta T
    """
}


