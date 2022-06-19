#! /bin/bash


## Accesing to working directory
cd /home/omicas/sarab/rnaseq_co_7/samples/bams

## Gene Expression Quantification
featureCounts -O -p -a ../../annotation/arabidopsis_thaliana.gtf -o counts.txt $(ls *.bam) -s 2 --verbose

