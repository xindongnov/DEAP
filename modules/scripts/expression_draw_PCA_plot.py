#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   expression_RNASeq_DE_analysis.py
@Time    :   2024/08/19 01:26:28
@Author  :   Xin Dong
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2024, XinDong
'''

import os
import sys
# import pickle as pkl
import pandas as pd
import numpy as np
import scipy 
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
from inmoose.pycombat import pycombat_norm, pycombat_seq
import argparse
from sklearn import decomposition
import matplotlib.pyplot as plt
import adjustText

plt.rcParams.update({
    'figure.figsize': [4, 4],
    'font.size' : 10,
    'font.family': 'Arial',
    'font.style' : 'normal',
    'font.weight':'normal',
    'figure.titleweight': 'normal',
    'axes.labelsize': 10 ,
    'axes.titleweight': 'normal',
    'axes.labelweight': 'normal',
    'axes.spines.right': False,
    'axes.spines.top': False,
})

def get_label_color(label):
    label_color = []
    n = 1
    for i,j in enumerate(label):
        if i == 0:
            label_color.append(n)
        else:
            if label[i] == label[i-1]:
                label_color.append(n)
            else:
                n += 1
                label_color.append(n)
    return label_color

def main():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser()
    parser.add_argument('-d', '--design', help='The design file', required=True)
    parser.add_argument('-c', '--count', help='the count path to save', default="count.txt")
    parser.add_argument('-o', '--condition_output', help='the output path', default='condition_PCA.pdf')
    parser.add_argument('-b', '--batch_output', help='the output path', default='batch_PCA.pdf')
    args = parser.parse_args()


    count_path = args.count
    design_path = args.design
    output_path = args.output
    
    design = pd.read_csv(design_path, sep='\t', index_col=0)
    # print(design)
    all_counts_df = pd.read_csv(count_path, sep='\t', index_col=0)
    # print(all_counts_df)

    ref_condition = design[design[compare] == 0].condition.unique()[0]
    treat_condition = design[design[compare] == 1].condition.unique()[0]
    # output_path = 'lcy/%s_vs_%s/lfc_%s_Pvalue_%s/' % (treat_condition, ref_condition, lfc_cutoff, pvalue_cutoff)
    # print(output_path)
    do_batch=True
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    # lfc_cutoff = 0.5
    # pvalue_cutoff = 1e-5

    print(ref_condition, treat_condition, output_path, lfc_cutoff, pvalue_cutoff)

    needed_sample = design.index[design[compare] != '']
    metadata = design.loc[needed_sample, ['condition', 'batch', compare]]
    # print(metadata)

    counts_df = all_counts_df.loc[needed_sample]
    counts_df = counts_df.drop(columns='-').dropna(axis=1)
    # print(counts_df)

    # # remove batch effect
    # for i in metadata['batch'].unique():
    #     if metadata[metadata['batch'] == i].__len__() == 1:
    #         do_batch=False
    # if do_batch:
    #     counts_df = pycombat_seq(counts_df.T,metadata['batch']).T

    # do PCA
    X_reduced = decomposition.PCA(n_components=2, svd_solver='arpack').fit_transform(counts_df)
    label_color = get_label_color(metadata['condition'])
    fig, ax= plt.subplots(figsize=(6, 6))
    scatter = ax.scatter(X_reduced[:, 0], X_reduced[:, 1], label=counts_df.index, c=label_color, s=40, cmap='tab20')
    texts = [plt.text(X_reduced[i, 0], X_reduced[i, 1], txt) for i, txt in enumerate(counts_df.index)]
    adjustText.adjust_text(texts, arrowprops=dict(arrowstyle='->', color='grey'))
    fig.savefig(os.path.join(output_path, 'condition_PCA.pdf'))

    fig, ax= plt.subplots(figsize=(6, 6))
    scatter = ax.scatter(X_reduced[:, 0], X_reduced[:, 1], label=counts_df.index, c=metadata['batch'].values, s=40, cmap='tab20')
    texts = [plt.text(X_reduced[i, 0], X_reduced[i, 1], txt) for i, txt in enumerate(counts_df.index)]
    adjustText.adjust_text(texts, arrowprops=dict(arrowstyle='->', color='grey'))
    fig.savefig(os.path.join(output_path, 'batch_PCA.pdf'))

if __name__ == '__main__':
    main()