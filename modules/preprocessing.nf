#!/usr/bin/env nextflow

nextflow.enable.dsl=2


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*                   TRIMMING adapter/quality trim                           */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process TRIMMING {
    publishDir "${params.outdir}/trimgalore", mode: 'copy', pattern: "*.html"

    input:
        tuple val(pair_id), path(reads)
        val trim_args

    output:
        tuple val(pair_id), path('*1.fq.gz'), path('*2.fq.gz'),     emit: trimread
        path "*fastqc.html" ,                                       emit: fastqc

    script:
    """
    fastqc -t $task.cpus ${reads[0]} ${reads[1]}
    trim_galore --cores $task.cpus --fastqc --gzip $trim_args --paired ${reads[0]} ${reads[1]}
    """
}


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*                    DEDUPE Remove duplicates                               */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

process DEDUPE {
    input:
        tuple val(pair_id), path(reads1), path(reads2)

    output:
        tuple val(pair_id), path('reads.unique.1.fq.gz'), path('reads.unique.2.fq.gz'),      emit: dedupreads

    script:
    """
    tally -i $reads1 -j $reads2 -o reads.unique.1.fq.gz -p reads.unique.2.fq.gz --pair-by-offset --with-quality
    """
}

///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*                    NORM normalize coverage                                */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////
process NORM {
    publishDir "${params.outdir}/normalized", mode: 'copy', pattern: "*.normalization_stats.txt" 

    input:
        tuple val(pair_id), path(reads1), path(reads2)

    output:
        tuple val(pair_id), path('reads.norm.1.fq.gz'), path('reads.norm.2.fq.gz'),         emit: normreads
        path "*.normalization_stats.txt",                                                   emit: stats

    script:
    """
    tadpole.sh overwrite=true mode=correct k=50 in=$reads1 in2=$reads2 out=reads.tadpole.1.fq.gz out2=reads.tadpole.2.fq.gz 
    bbnorm.sh overwrite=true mindepth=1 target=100 in=reads.tadpole.1.fq.gz in2=reads.tadpole.2.fq.gz out=reads.norm.1.fq.gz out2=reads.norm.2.fq.gz
    stats.sh overwrite=true in=reads.norm.1.fq.gz out=${pair_id}.normalization_stats.txt 
    """
}


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*                        CORRECT correct reads using kmers                  */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////
process CORRECT {

    input:
        tuple val(pair_id), path(reads1), path(reads2)

    output:
        tuple val(pair_id), path('reads.corrRep.1.fq.gz'), path('reads.corrRep.2.fq.gz'),         emit: correctreads

    script:
    """
    zcat $reads1 > reads.1.fq
    zcat $reads2 > reads.2.fq
    tadpole.sh in=reads.1.fq in2=reads.2.fq out=reads.1.corr.fq out2=reads.2.corr.fq mode=correct k=$params.correct_kmer ecc=t -Xmx1g prealloc=t prefilter=1 prepasses=auto
    repair.sh overwrite=true in1=reads.1.corr.fq in2=reads.2.corr.fq out1=reads.corrRep.1.fq.gz out2=reads.corrRep.2.fq.gz repair
    """
}

///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*                 MERGE merge forward and reverse reads overlapping         */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////
process MERGE {
    publishDir "${params.outdir}/merge", mode: 'copy', pattern: "*inserts.txt" 

    input:
        tuple val(pair_id), path(reads1), path(reads2)

    output:
        tuple val(pair_id), path('reads.merged.fq.gz'),                 emit: overlapreads
        path "*inserts.txt",                                            emit: inserts

    script:
    """
    bbmerge.sh in=$reads1 in2=$reads2 out=reads.merged.fq.gz outu1=reads1.unmerged.fq.gz outu2=reads2.unmerged.fq.gz outinsert=${pair_id}.inserts.txt
    """
}


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*                 SPLIT split based on tag at start of SE reads             */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////
process SPLIT {
    publishDir "${params.outdir}/split_tags", mode: 'copy', pattern: "*fq.gz" 

    input:
        tuple val(pair_id), path(reads)

    output:
        tuple val(pair_id), path('*.fq.gz'),                        emit: splitreads

    script:
    """
    $projectDir/bin/split_barcode.py --input $reads --prefix $pair_id --len $params.tag_len
    """
}



///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*                    CLUSTER cluster identical reads                        */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////
process CLUSTER {
    publishDir "${params.outdir}/clustering", mode: 'copy', pattern: "*fa" 

    input:
        tuple val(pair_id), path(reads)

    output:
        tuple val(pair_id), path('*.cluster_consensus.fa'),          emit: clusterreads

    script:
    """
    vsearch --threads $task.cpus --cluster_size $reads --strand both --consout ${pair_id}.cluster_consensus.fa --id 0.9
    """
}

