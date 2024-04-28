#!/bin/bash

# this bash file is used to generate rsem genome reference

# usage: bash reference_rsem_index.sh <reference_output_path> <fasta_path> <gtf_path>

REF_PATH=$1
FASTA_PATH=$2
GTF_PATH=$3

mkdir ${REF_PATH}
gunzip -c ${FASTA_PATH} > ${FASTA_PATH}.rsemtmp
gunzip -c ${GTF_PATH} > ${GTF_PATH}.rsemtmp

rsem-prepare-reference --gtf ${GTF_PATH}.rsemtmp ${FASTA_PATH}.rsemtmp ${REF_PATH}/rsem

rm ${FASTA_PATH}.rsemtmp
rm ${GTF_PATH}.rsemtmp