#!/bin/bash

zcat $1 | awk -v OFS="\t" '{print $1, $2, $3, $11}' - > $2
bedGraphToBigWig $2 /data/kzkarttu/annotations/genome/hg38.chrom.sizes $3

bgzip $2
tabix -p bed $2
