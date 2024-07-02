#!/bin/bash

/data/kzkarttu/software/modkit/modkit pileup  ../align/concatenated/GP5d_DNMTi_HDACi_concatenated.bam GP5d_DNMTi_HDACi_bedMethyl.bed \
    -t 32 \
    --log-filepath pileup.log \
    --ref /data/kzkarttu/annotations/genome/GRCh38.primary_assembly.genome.fa \
    --preset traditional

# Find all GCG motifs for filtering out ambiguous calls
# /data/kzkarttu/software/modkit/modkit motif-bed /data/kzkarttu/annotations/genome/GRCh38.primary_assembly.genome.fa  GCG 0 1> gcg_modifs.bed
