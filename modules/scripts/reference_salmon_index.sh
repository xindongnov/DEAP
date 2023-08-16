#!/bin/bash

# this bash file is used to generate STAR genome reference

# usage: bash reference_salmon_index.sh <reference_output_path> <fasta_path> <gtf_path> <threads>

REF_PATH=$1
FASTA_PATH=$2
GTF_PATH=$3
THREADS=$4

salmon index -t $FASTA_PATH --gencode $GTF_PATH -p $THREADS -i $REF_PATH