#!/bin/bash

#
# Extract reads in fastq format from bam file for PE reads where at least 1 read maps
#
# Inputs:
# $1: Input bam file prefix without file extension
# $2:Number of cores
# 
# Output: 
# $prefix.1.fq $prefix.2.fq, gzipped fastq files of reads that have at least 1 end mapping
#
function extract_bam_reads {
  local prefix=$1   
  local cores=$2
    
  if [[ $prefix == *".bam" ]]; then
    prefix=${prefix%.bam}
  fi
  samtools view --threads $cores -bh -F 4 -f 8 $prefix.bam > $prefix.1.bam
  samtools view --threads $cores -bh -F 8 -f 4 $prefix.bam > $prefix.2.bam
  samtools view --threads $cores -bh -F12 $prefix.bam > $prefix.3.bam
  samtools merge --threads $cores -f $prefix.merged.bam $prefix.1.bam $prefix.2.bam $prefix.3.bam 2> /dev/null
  samtools sort --threads $cores -n $prefix.merged.bam -o $prefix.merged.sorted.bam 2> /dev/null
  bedtools bamtofastq -i $prefix.merged.sorted.bam -fq $prefix.1.fq -fq2 $prefix.2.fq
  gzip *fq
}

#
# Extract reads in fastq format from bam file for SE reads
#
# Inputs:
# $1: Input bam file prefix without file extension
# $2:Number of cores
# 
# Output: 
# $prefix.fq, gzipped fastq file of mapped reads 
#
function extract_bam_reads_se {
  local prefix=$1   
  local cores=$2
    
  if [[ $prefix == *".bam" ]]; then
    prefix=${prefix%.bam}
  fi
  samtools view --threads $cores -bh -F 4 $prefix.bam > $prefix.map.bam
  #samtools sort --threads $cores -n $prefix.merged.bam -o $prefix.merged.sorted.bam 2> /dev/null
  bedtools bamtofastq -i $prefix.map.bam -fq $prefix.fq
  gzip *.fq
}

#
# Extract sequences from a fasta file based on a text file of positions of the
# format 'id' or 'id:start-stop' where id is the name of an element in the fasta file
#
# $1: file containing sequences to extract, position optional
# $2: genome file
# $3: output file
#
function extract_seq {
    local extract=$1
    local genome=$2
    local output=$3
    
    samtools faidx $genome
    >$output
    while IFS=' ' read col1 col2
    do
        if [ ! -z ${col2} ]; then 
            samtools faidx $genome $col1 | sed "s/>.*/>$col2/g" >> $output
        else 
            samtools faidx $genome $col1 >> $output
        fi
    done <$extract
}


#
# Extract blast hit sequences from a fasta file. 
# Scaffolds that blast on the negative strand are reverse complemented.
#
# $1: input file, blast output in format 6
# $2: sequence fasta file to extract from
# $3: output file
# $4: T/F, whether to reverse complement if hit is on the reverse strand
#
function extractBlastedScaff {
    local input=$1
    local sequence=$2
    local output=$3
    local stranded=$4
    
    if [ "$stranded" == "T" ]; then
        awk '$7 <= $8 {print $1}' $input |sort -k1,1 |uniq > $input.extractScaffForward.bed
        extract_seq $input.extractScaffForward.bed $sequence $output.forward.fa
        awk '$7 > $8 {print $1}' $input |sort -k1,1 |uniq > $input.extractScaffReverse.bed
        extract_seq $input.extractScaffReverse.bed $sequence $output.reverse.fa
        mv $output.forward.fa $output
        seqtk seq -r $output.reverse.fa >> $output
        rm $output.reverse.fa
    else
        awk '{print $1}' $input |sort -k1,1 |uniq > $input.bed
        extract_seq $input.bed $sequence $output
    fi
}

#
# Convert fasta file based on blast so that all seqeunces have the same strandedness.
#
function convert_strand {
    local sequences=$1             # fasta file of nucleotide sequences
    local strandfile=$2            # File containing strand info in format "name strand"
    local outfile=$3                # Name for output file

    >$outfile.forward.fa
    >$outfile.reverse.fa
    >$outfile
    samtools faidx $sequences
    while IFS=' ' read myseq strand
    do
        if [ $strand == "+" ]; then
            samtools faidx $sequences $myseq >> $outfile.forward.fa
        else
            samtools faidx $sequences $myseq >> $outfile.reverse.fa
        fi 
    done < $strandfile
    cut -d " " -f 1 $strandfile |sort -k1,1 |uniq  > $outfile.names
    if [ -s $outfile.reverse.fa ]; then
        seqtk seq -r $outfile.reverse.fa >> $outfile
    fi
    cat $outfile.forward.fa >>  $outfile 
}

#
# Blast the extracted sequence against the consensus to cluster based on similarity to them
#
function create_cluster {
    local consensus=$1            # Consensus sequences to compare to in fasta format
    local sequences=$2            # Sequences of interest in fasta format
    local dbtype=$3                # dbtype for blast, nucl or prot
    local outputdir=$4            # Output directory to write results to
    local cores=$5

    basecon=${consensus##*/}
    basecon=${basecon%.*}
    baseseq=${sequences##*/}
    baseseq=${baseseq%.*}

    mkdir -p $outputdir/${baseseq}_clusters
    samtools faidx $consensus

    cp $consensus $outputdir/$basecon.fa
    existing_files=$(ls $outputdir/${baseseq}_clusters/*consensus_clusters* 2> /dev/null |wc -l)
    if [ $existing_files -gt 0 ]; then
        rm $outputdir/${baseseq}_clusters/*consensus_clusters*
    fi
    fasta36 -E 1E-10 -T $cores -m 8 -b 1 $sequences $outputdir/$basecon.fa | awk '{if ($9 > $10 || $7 > $8){$11="-"} else{$11="+"}; print $0}' > $outputdir/$baseseq.consensus_blast_clustering.txt
    awk -v output=$outputdir/${baseseq}_clusters/ -v base=$baseseq '{print $1,$11 >> ""output"/"base".consensus_clusters."$2".txt" }' $outputdir/$baseseq.consensus_blast_clustering.txt

    for file in $outputdir//${baseseq}_clusters/$baseseq.consensus_clusters.*.txt
    do
        basename=${file##*/}
        conname=${basename##*consensus_clusters.}
        conname=${conname%.txt}
        cat $file |sort -k1,1 |uniq > $file.uniq
        convert_strand $sequences $file.uniq $file.stranded.fa
        >$outputdir/${baseseq}_clusters/$basename.fa
        extract_seq $file.stranded.fa.names $file.stranded.fa $outputdir/${baseseq}_clusters/$basename.fa
        samtools faidx $consensus $conname >> $outputdir/${baseseq}_clusters/$basename.fa
        muscle -clw -in $outputdir/${baseseq}_clusters/$basename.fa -out $outputdir/${baseseq}_clusters/$basename.muscle.clw \
            > $outputdir/${baseseq}_clusters/$basename.muscle.log 2> $outputdir/${baseseq}_clusters/$basename.muscle.err || echo "Error clustering $basename" 
    done
}

#
# getBRA
# Calculate best reciprocal aligments between fasta files, dna/dna, prot/prot, or prot/dna comparisons
#
# input: 
# $1: file1, 
# $2: file2, 
# $3: method, type of comparison, dna_dna prot_prot prot_dna
#
# output: 
# $file1.BRA.$method.txt: 
# $file1.BRA.full_list.$method.txt: 
# 
function getBRA {
  local file1=$1
  local file2=$2
  local method=$3

  if [[ "$method" == "prot_prot" || "$method" == "dna_dna" ]]; then
    # file order: query, library
    fasta36 -b 1 -d 0 -m 8 -E 0.001 -T 8 $file1 $file2 > $file1.annot_transcripts.txt
    fasta36 -b 1 -d 0 -m 8 -E 0.001 -T 8 $file2 $file1 > $file1.annot_transcripts2.txt
  elif [ "$method" == "prot_dna" ]; then
    fasty36 -b 1 -d 0 -m 8 -E 0.001 -T 8 $file1 $file2 > $file1.annot_transcripts.txt
    tfasty36 -b 1 -d 0 -m 8 -E 0.001 -T 8 $file2 $file1 > $file1.annot_transcripts2.txt
  fi
  awk '{print $1, $2}' $file1.annot_transcripts.txt |sort -k 1,1|uniq > $file1.annot_transcripts.txt.out
  awk '{print $2, $1}' $file1.annot_transcripts2.txt |sort -k 1,1|uniq > $file1.annot_transcripts2.txt.out
  join $file1.annot_transcripts.txt.out $file1.annot_transcripts2.txt.out | awk '{ if ($2==$3) print $1, $2 }' |sort -k 1,1 |uniq > $file1.BRA.$method.txt
  join <(grep ">" $file2 |sed 's/ .*//g' |sed 's/>//g' |sort -k 1,1 |uniq) <(awk '{print $2,$1}' $file1.BRA.$method.txt |sort -k 1,1 |uniq) -a 1 > $file1.BRA.full_list.$method.txt
}

function ortho_count {
  grep ">" $2 | awk '{print $1}' |sed 's/>//g' |sed 's/_.*//g'|sed 's/chr//g' |sort -k 1,1 |uniq > chrs.txt
  printf "Chr\tTotal\tWithOrhtologs\n" > $1.ortholog_counts_alignment.txt
  while read chr
  do
    printf "$chr\t" >> $1.ortholog_counts_alignment.txt
    myhits=$(grep "$chr" $2 |wc -l)
    printf "$myhits\t" >> $1.ortholog_counts_alignment.txt
    myhits=$(grep "$chr" $1 | wc -l)
    printf "$myhits\n" >> $1.ortholog_counts_alignment.txt
  done < chrs.txt
}
