#!/bin/bash

/data/kzkarttu/software/dorado/bin/dorado aligner /data/kzkarttu/annotations/genome/GRCh38.primary_assembly.genome.fa /data/Samples/nanopore_seq/GP5d_DNMti_HDACi/20240528_1026_3E_PAW72130_0a7d2900/bam_pass/ -o ./ -t 32 --emit-summary
