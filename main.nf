#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {read_fastq}                                    from './modules/filehandling'
include {TRIMMING; SKIPTRIM; DEDUPE; CORRECT; NORM}     from './modules/preprocessing'
include {NGMALIGN; COMPLEXITYFILTER; SPADESASSEM}       from './modules/modules'
include {EXTRACTBAM; BLASTFILTER; CYCLEASSEM}           from './modules/modules'
include {EXTRACTEXONS; FINDEXONS; CLUSTER}              from './modules/modules'

/*
========================================================================================
                         cycleassembler
========================================================================================
 #### Homepage / Documentation
 https://github.com/ntoda03/
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    Usage:

    The pipeline can be run as followed:

      nextflow run cycleassembler 

    Mandatory arguments:

    Input files
        --reads  [file]                 Directory containing input fastq files to analyze

    References
        --reference [file]              Homologous sequences of interest to focus on in fasta format
        --reference_type [str]          ['dna' or 'prot']
        --seeds [file]                  Seed contigs from a previous assembly of this data to use to 
                                        identify initial reads to use

    Optional Arguments:

    Directories
        -w                              Scratch working directory
        --outdir                        Diretory to store output files in

    Trimming
        --clip_r1 [int]                 Remove int bases from the start of paired end read 1 (default: 0)
        --clip_r1_end   [int]           Remove int bases from the end of paired end read 1 (default: 0)
        --clip_r2_start [int]           Remove int bases from the start of paired end read 2 (default: 0)
        --clip_r2_end   [int]           Remove int bases from the end of reverse paired end read 2 (default: 0)
        --skip_trimming [bool]          Skip the adapter trimming step (default: false)
        --skip_dedupe [bool]            Skip the read deduplication step  (default: false)

    Assembly
        --initial_scaffolds [file]      Seed scaffolds to use, otherwise will be based on mapping to reference

    Exon extraction
        --exons [file]                  Fasta file containing exon sequences. The exons will be mapped
                                        to the assembled contigs and the corresponding sequences will be
                                        extracted.

""".stripIndent()
}

///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*                             Parameter checks                              */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

// Show help message
params.help = false
if (params.help) {
    helpMessage()
    exit 0
}


///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*                                Workflow                                   */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

workflow {
    
    // Input reads are all paired fq and fq.gz files in input dir
    reads_ch = read_fastq(params.reads)
    ref_ch = file(params.reference)

    if( params.reference_type == 'nucl' ){
        fasta_command = "fasta36"} 
    else {
        fasta_command = "fastx36"}

    // Basic read processing
    def trim_args = ''
    if ( params.clip_r1 != 0 ){ trim_args += '--clip_R1 ' + params.clip_r1 }
    if ( params.three_prime_clip_r1 != 0 ){ trim_args += '--three_prime_clip_R1 ' + params.three_prime_clip_r1 }
    if ( params.clip_r2 != 0 ){ trim_args += '--clip_R2 ' + params.clip_r2 }
    if ( params.three_prime_clip_r2 != 0 ){ trim_args += '--three_prime_clip_R2 ' + params.three_prime_clip_r2 }
    if( ! params.skip_trimming ){
        TRIMMING(reads_ch, trim_args)
        trim_ch = TRIMMING.out.trimread
    }
    else{
        SKIPTRIM(reads_ch, trim_args)
        trim_ch = SKIPTRIM.out.skiptrim
    }
    if( ! params.skip_dedupe ){
        DEDUPE(trim_ch)
        dedupe_ch = DEDUPE.out.dedupreads
    }
    else{
        dedupe_ch = trim_ch}
    NORM(dedupe_ch)
    
    // Initial assembly to get seed sequences
    NGMALIGN(NORM.out.normreads, ref_ch)
    EXTRACTBAM(NGMALIGN.out.ngmbam)
    COMPLEXITYFILTER(EXTRACTBAM.out.extractread, '-f 0 -t 0 -u 0')
    SPADESASSEM(COMPLEXITYFILTER.out.filterread, '--cov-cutoff 1')
    seeds_ch = SPADESASSEM.out.assembly

    // Iteratively assemble reads 
    BLASTFILTER(seeds_ch,ref_ch,fasta_command)
    CYCLEASSEM(BLASTFILTER.out.filtercontigs,NORM.out.normreads,ref_ch,fasta_command,params.maxit )

    // optional extraction of exon sequences
    if( params.exons ){
        exons_ch = file(params.exons, checkIfExists: true)
        EXTRACTEXONS(CYCLEASSEM.out.cyclecontigs, exons_ch)
        FINDEXONS(EXTRACTEXONS.out.exonseqs, exons_ch)
        CLUSTER(EXTRACTEXONS.out.exonseqs, exons_ch)
    }
}