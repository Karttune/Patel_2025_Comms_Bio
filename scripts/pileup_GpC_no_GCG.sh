#!/bin/bash

in_file=../data/Nanopore/GP5d_DNMTi_HDACi/align/concatenated/GP5d_DNMTi_HDACi_concatenated.bam
out_file=../data/Nanopore/GP5d_DNMTi_HDACi/methylation/GP5d_DNMTi_HDACi_bedMethyl_GpC_no_GCG.bed.gz

echo $in_file
echo $out_file

/data/kzkarttu/software/modkit/modkit pileup $in_file stdout \
    -t 40 \
    --log-filepath pileup_GpC_no_GCG.log \
    --ref /data/kzkarttu/annotations/genome/GRCh38.primary_assembly.genome.fa \
    --ignore h \
    --motif GCT 1 \
    --motif GCA 1 \
    --motif GCC 1 | bgzip > $out_file

tabix -p bed $out_file
