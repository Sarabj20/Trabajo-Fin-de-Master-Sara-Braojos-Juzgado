#! /bin/bash

SMP=$1
ACC=$2

## Accesing to working directory
cd /home/omicas/sarab/rnaseq_co/samples/$SMP

## Sample quality control and read mapping to reference genome
fastqc ${ACC}.fastq.gz
STAR --runThreadN 4 --genomeDir ../../genome/ --readFilesIn $ACC.fastq.gz --outFileNamePrefix $SMP --readFilesCommand gunzip -c

## Generating sorted bam file
samtools sort -o $SMP.bam ${SMP}Aligned.out.sam
rm ${SMP}Aligned.out.sam
samtools index $SMP.bam
bamCoverage -bs 10 --normalizeUsing CPM --bam $SMP.bam -o $SMP.bw

## Transcript assembly
stringtie -G ../../annotation/arabidopsis_thaliana.gtf -o $SMP.gtf -l $SMP $SMP.bam

## Preparing merge list file for transcriptome merging
echo /home/omicas/sarab/rnaseq_co/samples/$SMP/$SMP.gtf >> ../../results/merge_list.txt

## Moving bam file to another folder in order to perform featureCounts
mv $SMP.bam ../bams
