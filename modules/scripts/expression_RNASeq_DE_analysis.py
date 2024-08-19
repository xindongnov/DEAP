#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   expression_RNASeq_DE_analysis.py
@Time    :   2024/08/19 01:26:28
@Author  :   Xin Dong
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2024, XinDong
'''

# import os
import sys
# import pickle as pkl
import pandas as pd
# import numpy as np
# import scipy 
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
# from pydeseq2.utils import load_example_data

import argparse
# from sklearn import decomposition
# import matplotlib.pyplot as plt
# import adjustText

def main():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser()
    parser.add_argument('-c', '--count', help='the raw count table path', default="count.txt")
    parser.add_argument('-d', '--design', help='The design file', required=True)
    parser.add_argument('-cmp', '--compare', help='the compare type', default='treat_vs_ref')
    parser.add_argument('-o', '--output', help='the output path', default='DEseq.txt')
    args = parser.parse_args()

    count_path = args.count
    design_path = args.design
    compare = args.compare
    output_path = args.output
    design = pd.read_csv(design_path, sep=',', index_col=0)
    # print(design)
    all_counts_df = pd.read_csv(count_path, sep='\t', index_col=0).T
    # print(all_counts_df)

    ref_condition = design[design[compare] == 0].condition.unique()[0]
    treat_condition = design[design[compare] == 1].condition.unique()[0]

    # output_path = 'lcy/%s_vs_%s/lfc_%s_Pvalue_%s/' % (treat_condition, ref_condition, lfc_cutoff, pvalue_cutoff)
    # print(output_path)
    # do_batch=True
    # if not os.path.exists(output_path):
    #     os.makedirs(output_path)
    # lfc_cutoff = 0.5
    # pvalue_cutoff = 1e-5

    # print(ref_condition, treat_condition, output_path, lfc_cutoff, pvalue_cutoff)

    needed_sample = design.index[design[compare] != '']
    metadata = design.loc[needed_sample, ['condition', 'batch', compare]]
    # print(metadata)

    counts_df = all_counts_df.loc[needed_sample]
    counts_df = counts_df.round().astype(int)
    # counts_df = counts_df.drop(columns='-').dropna(axis=1)
    # print(counts_df)

    genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
    counts_df = counts_df[genes_to_keep]

    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors="condition",
        ref_level=['condition',ref_condition],
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()    
    stat_res = DeseqStats(dds, contrast=['condition',treat_condition,ref_condition], inference=inference)
    stat_res.summary()
    # deseq_table_path = os.path.join(output_path)
    deseq_table = stat_res.results_df.sort_values(by='padj')
    deseq_table.to_csv(output_path, sep='\t')

    # up_gene = deseq_table[(deseq_table['log2FoldChange'] > lfc_cutoff) & (deseq_table['padj'] < pvalue_cutoff)].sort_values(by='log2FoldChange', ascending=False).index
    # down_gene = deseq_table[(deseq_table['log2FoldChange'] < -lfc_cutoff) & (deseq_table['padj'] < pvalue_cutoff)].sort_values(by='log2FoldChange').index

    # with open(os.path.join(output_path, 'up_gene.txt'), 'w+') as file:
    #     file.write('\n'.join(up_gene))
    # with open(os.path.join(output_path, 'down_gene.txt'), 'w+') as file:
    #     file.write('\n'.join(down_gene))

if __name__ == '__main__':
    main()