#!/bin/bash

/data/kzkarttu/software/modkit/modkit pileup  ../align/concatenated/GP5d_DNMTi_HDACi_concatenated.bam stdout \
    -t 40 \
    --log-filepath pileup_no_GCG.log \
    --ref /data/kzkarttu/annotations/genome/GRCh38.primary_assembly.genome.fa \
    --ignore h \
    --motif TCG 1 \
    --motif CCG 1 \
    --motif ACG 1 | bgzip > GP5d_DNMTi_HDACi_bedMethyl_no_GCG.bed.gz

tabix -p bed GP5d_DNMTi_HDACi_bedMethyl_no_GCG.bed.gz
