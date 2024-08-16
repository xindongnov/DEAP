#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   expression_RNASeq_get_matrix.py
@Time    :   2023/07/01 14:53:32
@Author  :   Xin Dong
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2023, XinDong
'''

import sys
import os
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
    # parser.add_argument('-t', '--type', help='type of data', choices = ['Rawcount', 'TPM'], required=True)
    args = parser.parse_args()

    input = args.input
    result_path = args.output
    ref = args.reference
    # dtype = args.type

    in_matrix = pd.read_csv(input, sep='\t', index_col=0)
    gtf_anno = pd.read_csv(ref, sep='\t')
    index = list(set(in_matrix.index).intersection(gtf_anno['transcript_id']))
    in_matrix = in_matrix.loc[index,:].copy()
    # gtf_anno = gtf_anno[gtf_anno['transcript_id'].isin(index)].copy()
    mapping_dict = gtf_anno.set_index('transcript_id').loc[in_matrix.index, 'gene_name'].to_dict()
    in_matrix['name'] = in_matrix.index.map(mapping_dict)

    output = in_matrix.groupby('name').max()
    output.to_csv(result_path, sep='\t', index=True, header=True)

if __name__ == "__main__":
    main()