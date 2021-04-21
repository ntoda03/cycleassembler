
# poleanalyse/cycleassembler

**Iterative assembly pipeline that extends from generated seed contigs**.

## Introduction

The **poleanalyse/cycleassembler** pipeline takes as input illumina reads and reference sequences of interest. The pipeline assembles initial contigs that map to the reference and then uses these new contigs to itertively search for more reads to extend the contigs. This pipeline is best used for plastid genome assembly, low coverage genome skimming assembly of histone genes, or assembly of target genes from whole genome data. The assembly method used is not memory efficient so large references files (e.g. whole eukaryote genomes) should not be used.

The pipeline uses [Nextflow](https://www.nextflow.io) to facilitate parallel and reproducible analyses. It uses Docker to make it easier to use and independent of software requirements.

### Requirements

This pipeline can be run with only [Docker](https://docs.docker.com/engine/installation/) and [Nextflow](https://www.nextflow.io) and an internet connection. A docker container containing all other dependencies will be downloaded by the pipeline and can be found [here](https://hub.docker.com/repository/docker/poleanalyse/cycleassembler). If Docker is not used then the individual dependencies must be installed.

* Required software
    - [Nextflow](https://www.nextflow.io)
* Using with Docker. Highly recommended.
    - [Docker](https://docs.docker.com/engine/installation/) must be installed 
    - Docker will handle all other dependencies when the pipeline is run with the option "-profile docker"
* Using without Docker. The following software must be installed and in the path:
    - python3 with pandas
    - trim-galore
    - bbmap
    - vsearch
    - nextgenmap
    - samtools
    - fastp
    - spades
    - seqtk
    - muscle
    - reaper
    - bedtools
    - fasta

### Quick start

Simply pull the project from github.

`nextflow pull ntoda03/cycleassembler -r master`

Then it can be run passing in only the reads to analyze and the reference genome.

`nextflow run cycleassembler -profile docker --reads "path/to/reads/*_R{1,2}.fq.gz" --reference path/to/sequences/sequences.fa`

### Citing this pipeline

This pipeline is part of the work of the pole analyse UMS2700 at the National Museum of Natural History in Paris.

If you use this pipeline please cite it using the DOI 10.5281/zenodo.4609123

Thanks!

## Running the pipeline

### Pipeline Summary

By default, the pipeline does the following:

Basic read processing
* Read trimming for adapters and low quality bases (`trimgalore`)
* Deduplicate identical reads (`tally`)
* Normalize high coverage regions (`bbnorm`)

Generate seed sequences
 * Align reads to reference with ~30% max divergence (`NGM`)
 * Extract aligned reads (`samtools`)
 * Filter low-complexity sequences (`fastp`)
 * Assemble extracted reads (`spades`)
 * Filter for seed sequences that align to reference (`fasta36`)

Iterative assembly
 * Align reads to seed sequences (`NGM`)
 * Extract aligned reads (`samtools`)
 * Filter low-complexity sequences (`fastp`)
 * Assemble extracted reads (`spades`)
 * Filter for contigs that map to reference (`fasta36`)
 * Repeat with the assembled contigs as the new seed sequences (default 5 cycles)

[optional] Extract precise reference matches from the final contigs. For use if using exon sequences as the reference sequences. Using gene sequences will cause issues at intron boundaries
 * Get contigs in correct orientation relative to reference (`seqtk`,`fastq36`)
 * Filter for alignments that span 80% of the exon length
 * Extract the sequences (`samtools`)
 * Get the best reciprical alignment between exons and sequences to have a single hit per reference exon (`fasta36`)
 * Cluster and align the hits (`samtools`,`fasta36`,`muscle`)

### How to use

```
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
    --outdir                        Directory to store results files in (default: ./results)

Pipeline control
    -resume                         Continue a previously running analysis that did not finish.
    
Profiles
    -profile docker                 Use docker to handle all the software dependencies. 

Reads
    --single_end [bool]             Whether reads are single end (true/false, default: false)
                                    This is not recommended since the benefit of this pipeline comes
                                    mostly from finding paired reads where only one end maps to extend
                                    the contigs.

Trimming options
    --clip_r1 [int]                 Remove int bases from the start of paired end read 1 (default: 0)
    --clip_r1_end   [int]           Remove int bases from the end of paired end read 1 (default: 0)
    --clip_r2_start [int]           Remove int bases from the start of paired end read 2 (default: 0)
    --clip_r2_end   [int]           Remove int bases from the end of reverse paired end read 2 (default: 0)
    --skip_trimming [bool]          Skip the adapter trimming step (true/false, default: false)
    --skip_dedupe [bool]            Skip the read deduplication step (true/false, default: false)
    --output_trimmed [bool]         Ouput the trimmed reads to the results directory (true/false, default: false)

Assembly 
    --maxit [int]                   The number of cycles to run for the iterative assembly. Increase this if
                                    you need to travers divergent intergenic or intronic sequences between 
                                    genes or exons. (default: 5)
    --orient [bool]                 Whether to orient contigs relative to the reference (true/false, default: true)
                                    This should be turned off if a contig may have matches in multiple orientations

Exon extraction
    --exons [file]                  Fasta file containing exon sequences. The exons will be mapped
                                    to the assembled contigs and the corresponding sequences will be
                                    extracted.
```

### Output files

In the output folder the following output folders and files will be created.

* FASTQC_raw_reads/ 
    - Quality control information of raw reads before trimming
* FASTQC_trimmed_reads/ 
    - Quality control information of reads after trimming
* trimmed_reads/ [optional] 
    - Trimmed reads in gzipped fastq format
* normalized
    - Some basic stats on reads, generally not interesting
* assembled_contigs
    - A fasta file containing the sequences assembled by the program
* exons
    - Optional sequences extracted that match the provided exon sequences
* exon_clusters
    - Each provided exon is matched to a best reciprocal hit of a sequences assembled if possible
