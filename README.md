
# poleanalyse/cycleassembler

**Iterative assembly pipeline that extends from generated seed contigs**.

## Introduction

The **poleanalyse/cycleassembler** pipeline takes as input illumina reads and reference sequences of interest. The pipeline assembles initial contigs that map to the reference and then uses these new contigs to itertively search for more reads to extend the contigs. This pipeline is best used for plastid genome assembly, low coverage genome skimming assembly of histone genes, or assembly of target genes from whole genome data. The assembly method used is not memory efficient so large references files (e.g. whole eukaryote genomes) should not be used.

The pipeline uses [Nextflow](https://www.nextflow.io) to facilitate parallel and reproducible analyses. It uses Docker to make it easier to use and independent of software requirements.

### Requirements

This pipeline can be run with only [Docker](https://docs.docker.com/engine/installation/) and [Nextflow](https://www.nextflow.io) and an internet connection. A docker container containing all other dependencies will be downloaded by the pipeline and can be found [here](https://hub.docker.com/repository/docker/poleanalyse/cycleassembler).

### Quick start

Simply pull the project from github.

`nextflow pull ntoda03/cycleassembler -r master`

Then it can be run passing in only the reads to analyze and the reference genome.

`nextflow run cycleassembler -profile docker --reads "*_R{1,2}.fq.gz" --reference sequences.fa`

### Citing this pipeline

If you use this pipeline please cite it using the DOI 10.5281/zenodo.4609123

## Running the pipeline

### Pipeline Summary

By default, the pipeline does the following:

Basic read processing
* Read trimming (`trimgalore`)
* Deduplicate reads (`tally`)
* Normalize coverage (`bbnorm`)

Generate seed sequences
 * Align reads to reference (`NGM`)
 * Extract aligned reads (`samtools`)
 * Filter low-complexity sequences (`fastp`)
 * Assemble reads (`spades`)
 * Filter for seed sequences that map to reference (`fasta36`)

Iterative assembly
 * Align reads to seed sequences (`NGM`)
 * Extract aligned reads (`samtools`)
 * Filter low-complexity sequences (`fastp`)
 * Assemble reads (`spades`)
 * Filter for contigs that map to reference (`fasta36`)
 * Repeat with the assembled contigs as the new seed sequences

[optional] Extract precise reference matches from the final contigs. For use if using exon sequences as the reference sequences. Using gene sequences will cause issues at intron boundaries
 * Get contigs in correct orientation (`seqtk`,`fastq36`)
 * Filter for hits that span 80% of the exon length
 * Extract the sequences (`samtools`)
 * Get the best reciprical alignment between exons and sequences to have a single hit per reference exon (`fasta36`)
 * Cluster and align the hits (`samtools`,`fasta36`,`muscle`)

### How to use

```
The pipeline can be run as followed:

  nextflow run cycleassembler [options] --reads "reads_{1,2}.fq.gz" --reference sequences.fa

Mandatory arguments:

Input files
    --reads  [file]                 Input paired end fastq files to analyze in quotes. 
                                    {1,2} indicates the read pairing number and this must  come after the id.
                                    A semicolon separated list and wildcards are be accepted.

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

Exon extraction
    --exons [file]                  Fasta file containing exon sequences. The exons will be mapped
                                    to the assembled contigs and the corresponding sequences will be
                                    extracted.
```

### Output files

In the output folder the following output folders and files will be created.

* trimgalore/ 
    - Quality control information on reads
* normalized
    - Some basic stats on reads
* assembled_contigs
    - A fasta file containing the sequences assembled by the program
* exons
    - Optional sequences extracted that match the provided exon sequences
* exon_clusters
    - Each provided exon is matched to a best reciprocal hit of a sequences assembled if possible
