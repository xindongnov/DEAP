#!/bin/bash

# this bash file is used to generate STAR genome reference

# usage: bash reference_STAR_index.sh <reference_output_path> <fasta_path> <gtf_path> <threads>

REF_PATH=$1
FASTA_PATH=$2
GTF_PATH=$3
THREADS=$4

gunzip -c $FASTA_PATH > ${FASTA_PATH}.tmp
gunzip -c $GTF_PATH > ${GTF_PATH}.tmp

STAR --runThreadN $THREADS \
--runMode genomeGenerate \
--genomeDir $REF_PATH \
--genomeFastaFiles ${FASTA_PATH} \
--sjdbGTFfile ${GTF_PATH} \
--limitGenomeGenerateRAM 104426985056 \
--sjdbOverhang 100

rm ${FASTA_PATH}.tmp
rm ${GTF_PATH}.tmp