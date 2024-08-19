#!/bin/bash

# this bash file is used to generate salmon genome reference

# usage: bash reference_salmon_index.sh <reference_output_path> <fasta_path> <threads>

REF_PATH=$1
FASTA_PATH=$2
THREADS=$3

salmon index -t ${FASTA_PATH} -p ${THREADS} -i ${REF_PATH}.nosplit

salmon index -t ${FASTA_PATH} -p ${THREADS} -i ${REF_PATH} --gencode