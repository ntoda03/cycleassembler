#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {TRIMMING; DEDUPE; NORM}                    from './modules/preprocessing'
include {SEEDASSEMBLY; CYCLEASSEMBLER}              from './modules/cycleassembly'
include {read_pefastq_dir}              from './modules/filehandling'

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
        --input_dir  [dir]              Directory containing input fastq files to analyze
        --dir_format [str]              Format of files in the input directory [pe', 'interleaved']
        --pe1 [file/list]               File or comma separated list of input paired end fastq files, read 1
        --pe2 [file/list]               File or comma separated list of input paired end fastq files, read 2
        --interleaved [file/list]       File or comma separated list of input interleaved fastq files

    References
        --reference [file]              Homologous sequences of interest to focus on in fasta format
        --reference_type [str]          ['dna' or 'prot']

    Optional Arguments:

    Directories
        -w                              Scratch working directory
        --outdir                        Diretory to store output files in

    Trimming
        --clip_r1_start [int]           Remove int bases from the start of paired end read 1 (default: 0)
        --clip_r1_end   [int]           Remove int bases from the end of paired end read 1 (default: 0)
        --clip_r2_start [int]           Remove int bases from the start of paired end read 2 (default: 0)
        --clip_r2_end   [int]           Remove int bases from the end of reverse paired end read 2 (default: 0)
        --skip_trimming [bool]          Skip the adapter trimming step (default: false)

    Assembly
        --initial_scaffolds [file]      Seed scaffolds to use, otherwise will be based on mapping to reference

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

////////////////////////////////////////////////////
/*                Default parameters              */
////////////////////////////////////////////////////

// files
params.outdir = './'
params.input_dir = '/home/nick/Bureau/Programming/nextflow/data/ggal/'

// trimming
params.clip_r1_start = 0
params.clip_r1_end = 0
params.clip_r2_start = 0
params.clip_r2_end = 0
params.skip_trimming = false

// cycle assmebly
params.reference = "/home/nick/Bureau/Programming/nextflow/data/ggal/transcriptome.fa"
params.reference_type = "nucl"
params.initial_scaffolds = ''
params.maxit = 10


////////////////////////////////////////////////////
/*                 Validate inputs                */
////////////////////////////////////////////////////

if( params.pe1 ){ params.pe1 = params.pe1.split(',') }
if( params.pe2 ){ params.pe2 = params.pe2.split(',')}

pathParamList = [ params.reference, params.pe1, params.pe2 ].collect()
for( param in pathParamList ){ if(param){ file(param,checkIfExists:true) } }

///////////////////////////////////////////////////////////////////////////////
/*                                                                           */
/*                                Workflow                                   */
/*                                                                           */
///////////////////////////////////////////////////////////////////////////////

workflow {
    
    //params.outdir += '/cycleassembler/'

    // Input reads are all paired fq and fq.gz files in input dir
    reads_ch = read_pefastq_dir(params.input_dir)
    
    def trim_args = ''
    if ( params.clip_r1_start != 0 ){ trim_args += '--clip_R1 ' + params.clip_r1_start }
    if ( params.clip_r1_end != 0 ){ trim_args += '--three_prime_clip_R1 ' + params.clip_r1_end }
    if ( params.clip_r2_start != 0 ){ trim_args += '--clip_R2 ' + params.clip_r2_start }
    if ( params.clip_r2_end != 0 ){ trim_args += '--three_prime_clip_R2 ' + params.clip_r2_end }

    TRIMMING(reads_ch, trim_args)
    DEDUPE(TRIMMING.out.trimread)
    NORM(DEDUPE.out.dedupreads)
    SEEDASSEMBLY(NORM.out.normreads, params.reference)
    CYCLEASSEMBLER(SEEDASSEMBLY.out, params.reference, params.reference_type, params.maxit)
    //CYCLEASSEMBLER.out.view()
}
