#!/bin/bash

# this bash file is used to generate tx2genename and tx2geneid from gtf
# defalut the column of tx_id, gene_id, gene_name is 12, 10, 16
# usage: bash reference_tx2genemap.sh <gtf_file>

GTF_file=$1
# tx_id_col=$2
# gene_id_col=$3
# gene_name_col=$4


zcat $GTF_file | awk '$3 == "transcript" {print $12"\t"$16}' | tr -d '";' > tx2genename.tsv
zcat $GTF_file | awk '$3 == "transcript" {print $12"\t"$10}' | tr -d '";' > tx2geneid.tsv