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
import argparse
import pandas as pd
import numpy as np

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

def volcano_plot(deseq_table, lfc_cutoff=0.5, pvalue_cutoff=1e-5, figsize=(5,5), add_text='auto', n_gene_name=15, plot_gene_names = [],
                 fig_title='', save_path='volcano_plot.pdf'):
    '''
    add_text: 'auto' or 'define' or 'none'
    '''
    result = pd.DataFrame()
    result['x'] = deseq_table['log2FoldChange']
    result['y'] = -np.log10(deseq_table['padj'])
    
    # set cutoff
    x_threshold = lfc_cutoff
    y_threshold = -np.log10(pvalue_cutoff)
    
    # group to up, normal, down
    result['group'] = 'dimgrey'
    result.loc[(result.x > x_threshold)&(result.y > y_threshold),'group'] = 'tab:red' # up regulated
    result.loc[(result.x < -x_threshold)&(result.y > y_threshold),'group'] = 'tab:blue' # down regulated
    # set which gene use to plot gene name
    if add_text == 'auto':
        result['xy'] = abs(result['x'] * result['y']) 
        plot_gene_names = result[result['group'] != 'dimgrey'].sort_values('xy', ascending=False).head(n_gene_name).index
    
    # plot
    fig, ax = plt.subplots(figsize=figsize) 
    ax.set(title=fig_title)
    ax.scatter(result['x'], result['y'], s=2, c=result['group'])
    ax.set_ylabel('-Log10(P adj.)',fontweight='bold', fontsize=10)
    ax.set_xlabel('Log2 (fold change)',fontweight='bold', fontsize=10)
    ax.spines['right'].set_visible(False) 
    ax.spines['top'].set_visible(False) 

    plt.axvline(-x_threshold, color='dimgrey',linestyle='dashed', linewidth=1) 
    plt.axvline(x_threshold, color='dimgrey',linestyle='dashed', linewidth=1) 
    plt.axhline(y_threshold, color='dimgrey',linestyle='dashed', linewidth=1) 
    if add_text == 'auto' or add_text == 'define':
        texts = [plt.text(result.loc[txt, 'x'], result.loc[txt, 'y'], txt) for i, txt in enumerate(plot_gene_names)]
        adjustText.adjust_text(texts, arrowprops=dict(arrowstyle='->', color='grey'))
    fig.tight_layout()
    fig.savefig(save_path)


def main():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser()
    parser.add_argument('-i', '--input', help='the DEseq table path', default="DE.txt")
    parser.add_argument('-o', '--output', help='the output path', default='volcano_plot.pdf')
    parser.add_argument('-l', '--lfc', help='the lfc cutoff', default=0.5)
    parser.add_argument('-f', '--fdr', help='the fdr cutoff', default=1e-5)
    parser.add_argument('-a', '--add_genename', help='add gene name', choices=['auto','define','none'], default='auto')
    parser.add_argument('-n', '--n_gene', help='the n gene name', default=15)
    parser.add_argument('-g', '--gene_names', help='gene names, sep by comma', default='')
    parser.add_argument('-t', '--title', help='the fig title', default='')

    args = parser.parse_args()

    deseq_table = pd.read_csv(args.input, sep='\t', index_col=0)
    lfc_cutoff = float(args.lfc)
    pvalue_cutoff = float(args.fdr)
    add_text = args.add_genename
    n_gene_name = int(args.n_gene)
    plot_gene_names = args.gene_names.split(',')
    fig_title = args.title
    save_path = args.output

    volcano_plot(deseq_table, 
                 lfc_cutoff=lfc_cutoff, 
                 pvalue_cutoff=pvalue_cutoff, 
                 figsize=(4,4), 
                 add_text=add_text, 
                 n_gene_name=n_gene_name, 
                 plot_gene_names = plot_gene_names,
                 fig_title=fig_title, 
                 save_path=save_path)

if __name__ == '__main__':
    main()