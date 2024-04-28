#!/bin/bash

# this bash file is used to generate STAR genome reference

# usage: bash reference_STAR_index.sh <reference_output_path> <fasta_path> <gtf_path> <threads>

REF_PATH=$1
FASTA_PATH=$2
GTF_PATH=$3
THREADS=$4

gunzip -c ${FASTA_PATH} > ${FASTA_PATH}.startmp
gunzip -c ${GTF_PATH} > ${GTF_PATH}.startmp

STAR --runThreadN $THREADS \
--runMode genomeGenerate \
--genomeDir $REF_PATH \
--genomeFastaFiles ${FASTA_PATH}.startmp \
--sjdbGTFfile ${GTF_PATH}.startmp \
--limitGenomeGenerateRAM 104426985056 \
--sjdbOverhang 100

rm ${FASTA_PATH}.startmp
rm ${GTF_PATH}.startmp