
# poleanalyse/cycleassembler

## UNDER DEVELOPMENT

**Iterative assembly pipeline that extends from generated seed contigs**.

## Introduction

The **poleanalyse/cycleassembler** pipeline takes as input illumina reads and reference sequnces of interest. The pipeline assembles initial contigs that map to the reference and then uses these new contigs to itertively search for more reads to extend the contigs. 

The pipeline uses [Nextflow](https://www.nextflow.io) to facilitate parallel and reproducible analyses. It uses Docker to make it easier to use and independent of software requirements.

## Requirements

This pipeline can be run with only [Docker](https://docs.docker.com/engine/installation/) and [Nextflow](https://www.nextflow.io) and an internet connection. A docker container containing all other dependencies will be downloaded by the pipeline and can be found [here](https://hub.docker.com/repository/docker/poleanalyse/cycleassembler).

## Quick start


## Pipeline Summary

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

Extract reference matches from the final contigs [optional]
For use if using exon sequences as the reference sequences
Using gene sequences will cause issues at intron boundaries
 * Get contigs in correct orientation (`seqtk`,`fastq36`)
 * Filter for hits that span 80% of the exon length
 * Extract the sequences (`samtools`)
 * Get the best reciprical alignment between exons and sequences to have a single hit per reference exon (`fasta36`)
 * Cluster and align the hits (`samtools`,`fasta36`,`muscle`)

