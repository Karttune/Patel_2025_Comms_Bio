#!/bin/bash

# Konsta Karttunen, November 2024.
# Last update: 8.11.2024.

# A script for basic NaNOMe-seq (R10) data processing with Dorado and Modkit. This script aligns the reads from unaligned Nanopore .bam file output to the reference genome with Dorado and extracts the CG and GC methylation tables from the .bam files directly with Modkit.

# Sofware requirements:
# Dorado
# Modkit
# Samtools
# Tabix

# NB! The last steps for creating .bedgraph and .bw files require the accompanying mtsv2bedGraph.py script!
# The script can be downloaded from https://github.com/isaclee/nanopore-methylation-utilities

MTSV2BEDGRAPH_PATH=/data/kzkarttu/nanopore/nanopore-methylation-utilities/mtsv2bedGraph.py

# User arguments:
# 1. Directory of the unaligned .bam output files from Nanopore (/[RUN ID]/bam_pass/)
# 2. Output directory
# 3. Basename
# 4. Threads

BAM_INPUT_DIR=$1
OUTPUT_DIR=$2
BASENAME=$3
THREADS=$4

##### The filepaths below for the annotation files and software have to be entered:

# Reference genome .fasta file:
REF_PATH=/data/kzkarttu/annotations/genome/GRCh38.primary_assembly.genome.fa
# PATH for the Dorado executable:
DORADO_PATH=/data/kzkarttu/software/dorado/bin/dorado
# PATH for the Modkit executable:
MODKIT_PATH=/data/kzkarttu/software/modkit/modkit

set -ue

echo -e "\nStarting alignment...\n"

mkdir -p ${OUTPUT_DIR}/bams
mkdir -p ${OUTPUT_DIR}/merged
mkdir -p ${OUTPUT_DIR}/logs
mkdir -p ${OUTPUT_DIR}/methylation
mkdir -p ${OUTPUT_DIR}/bedgraphs

$DORADO_PATH aligner $BAM_INPUT_DIR -o ${OUTPUT_DIR}/bams -t $THREADS --emit-summary

echo -e "\nMerging alignment files...\n"
samtools merge ${OUTPUT_DIR}/merged/${BASENAME}_concatenated.bam ${OUTPUT_DIR}/bams/*.bam

echo -e "\nExtracting CpG methylation...\n"

$MODKIT_PATH pileup  ${OUTPUT_DIR}/merged/${BASENAME}_concatenated.bam stdout \
    -t $THREADS \
    --log-filepath ${OUTPUT_DIR}/logs/${BASENAME}_pileup_CpG.log \
    --ref $REF_PATH \
    --ignore h \
    --motif TCG 1 \
    --motif CCG 1 \
    --motif ACG 1 | bgzip > ${OUTPUT_DIR}/methylation/${BASENAME}_CpG_bedMethyl.bed.gz

tabix -p bed ${OUTPUT_DIR}/methylation/${BASENAME}_CpG_bedMethyl.bed.gz

echo -e "\nExtracting GpC methylation...\n"

$MODKIT_PATH pileup ${OUTPUT_DIR}/merged/${BASENAME}_concatenated.bam stdout \
    -t $THREADS \
    --log-filepath ${OUTPUT_DIR}/logs/${BASENAME}_pileup_GpC.log \
    --ref $REF_PATH \
    --ignore h \
    --motif GCT 1 \
    --motif GCA 1 \
    --motif GCC 1 | bgzip > ${OUTPUT_DIR}/methylation/${BASENAME}_GpC_bedMethyl.bed.gz

tabix -p bed ${OUTPUT_DIR}/methylation/${BASENAME}_GpC_bedMethyl.bed.gz

zcat ${OUTPUT_DIR}/methylation/${BASENAME}_CpG_bedMethyl.bed.gz | awk -v OFS="\t" '{print $1, $2, $3, $11}' - > ${OUTPUT_DIR}/bedgraphs/${BASENAME}_CpG_bedMethyl.bedgraph
zcat ${OUTPUT_DIR}/methylation/${BASENAME}_GpC_bedMethyl.bed.gz | awk -v OFS="\t" '{print $1, $2, $3, $11}' - > ${OUTPUT_DIR}/methylation/${BASENAME}_GpC_bedMethyl.bedgraph

bedGraphToBigWig ${OUTPUT_DIR}/bedgraphs/${BASENAME}_CpG_bedMethyl.bedgraph $CHROM_SIZES ${OUTPUT_DIR}/bedgraphs/${BASENAME}_CpG_bedMethyl.bw
bedGraphToBigWig ${OUTPUT_DIR}/bedgraphs/${BASENAME}_GpC_bedMethyl.bedgraph $CHROM_SIZES ${OUTPUT_DIR}/bedgraphs/${BASENAME}_GpC_bedMethyl.bw

bgzip ${OUTPUT_DIR}/bedgraphs/${BASENAME}_CpG_bedMethyl.bedgraph
bgzip ${OUTPUT_DIR}/bedgraphs/${BASENAME}_GpC_bedMethyl.bedgraph

tabix -p bed ${OUTPUT_DIR}/bedgraphs/${BASENAME}_CpG_bedMethyl.bedgraph
tabix -p bed ${OUTPUT_DIR}/bedgraphs/${BASENAME}_GpC_bedMethyl.bedgraph

