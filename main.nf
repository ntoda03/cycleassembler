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
    poleanalyse/cycleassembler version 1.0.1

    Usage:
    
The pipeline can be run as followed:

  nextflow run cycleassembler [options] --reads "path/to/reads/reads_{1,2}.fq.gz" --reference path/to/sequences/sequences.fa

Primary arguments:

Input files
    --reads  [file]                 Input paired end fastq files to analyze in quotes. Only paired data is accepted.
                                    {1,2} indicates the read pairing number and this must come after the id.
                                    A semicolon separated list and wildcards are be accepted, for eaxmple
                                    "*{1,2}.fq.gz" indicates all the paired fastq files in the directory and
                                    "reads1_{1,2}.fq.gz;reads2_{1,2}.fq.gz;reads3_{1,2}.fq.gz" is for 3 sets of paired files.
                                    This should include the path to the files if they are not in the current directory.

References
    --reference [file]              Homologous sequences of interest to focus on in fasta format.
                                    This should include the path to the files if they are not in the current directory.
    --reference_type [str]          Type fo reference sequence ['nucl' or 'prot', default: 'nucl')

Optional Arguments:

Directories
    -w                              Scratch working directory for temporary files (default: ./work)
    --outdir                        Diretory to store results files in (default: ./results)

Pipeline control
    -resume                         Continue a previously running analysis that did not finish.

Profiles
    -profile docker                 Use docker to handle all the software dependencies. 

Trimming options
    --clip_r1 [int]                 Remove int bases from the start of paired end read 1 (default: 0)
    --clip_r1_end   [int]           Remove int bases from the end of paired end read 1 (default: 0)
    --clip_r2_start [int]           Remove int bases from the start of paired end read 2 (default: 0)
    --clip_r2_end   [int]           Remove int bases from the end of reverse paired end read 2 (default: 0)
    --skip_trimming [bool]          Skip the adapter trimming step (default: false)
    --skip_dedupe [bool]            Skip the read deduplication step  (default: false)
    --output_trimmed [bool]         Ouput the trimmed reads to the results directory

Assembly seed sequences
    --seeds [file]                  Seed contigs from a previous assembly of this data to use to identify initial
                                    reads to use instead of doing an initial assembly. This is only necessary
                                    if no suitable reference is available.
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
    // input fasta reference also
    ref_ch = file(params.reference, checkIfExists: true)

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
        trim_ch = reads_ch
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
