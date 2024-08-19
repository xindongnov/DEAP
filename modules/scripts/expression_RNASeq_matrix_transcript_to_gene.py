#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   expression_RNASeq_get_matrix.py
@Time    :   2023/07/01 14:53:32
@Author  :   Xin Dong
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2023, XinDong
'''

# this script is a very simple script to transform the transcript id to gene id
# usage: python expression_RNASeq_matrix_transcript_to_gene.py -i <input matrix file> -o <output file> -r <reference gtf file> -f <from type> -t <to type>
# I only use the maximum expression value of the same transcript id as the gene id
# beacause the salmon and RSEM both give the gene estimate value, so this script simply not be called
# recommand to use the tximport or pytximport to perform the transformation

import sys
# import os
import pandas as pd
import argparse

def main():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser()
    parser.add_argument('-i', '--input', help='input matrix file', required=True)
    parser.add_argument('-o', '--output', help='the path to save', default="./output.txt")
    parser.add_argument('-r', '--reference', help='transformed gtf reference', required=True)
    parser.add_argument('-f', '--fromtype', help='to format', default="transcript_id", choices = ['transcript_id', 'gene_id'], required=True)
    parser.add_argument('-t', '--totype', help='to format', default="gene_name", choices = ['gene_name', 'gene_id'], required=True)
    args = parser.parse_args()

    input = args.input
    result_path = args.output
    ref = args.reference
    tfrom = args.fromtype
    tto = args.totype

    in_matrix = pd.read_csv(input, sep='\t', index_col=0)
    gtf_anno = pd.read_csv(ref, sep='\t')
    index = list(set(in_matrix.index).intersection(gtf_anno['transcript_id']))
    in_matrix = in_matrix.loc[index,:].copy()
    mapping_dict = gtf_anno.set_index(tfrom).loc[in_matrix.index, tto].to_dict()
    in_matrix['name'] = in_matrix.index.map(mapping_dict)

    output = in_matrix.groupby('name').max()
    output.to_csv(result_path, sep='\t', index=True, header=True)

if __name__ == "__main__":
    main()