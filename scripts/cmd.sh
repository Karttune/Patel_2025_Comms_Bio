#!/bin/bash

/data/kzkarttu/software/dorado/bin/dorado aligner /data/kzkarttu/annotations/genome/GRCh38.primary_assembly.genome.fa /data/Samples/nanopore_seq/OE19_DNMTi_HDACi/20240528_1026_3G_PAW73334_c0516cc8/bam_pass/ -o ./ -t 32 --emit-summary
