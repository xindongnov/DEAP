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

# def get_label_color(label):
#     label_color = []
#     n = 1
#     for i,j in enumerate(label):
#         if i == 0:
#             label_color.append(n)
#         else:
#             if label[i] == label[i-1]:
#                 label_color.append(n)
#             else:
#                 n += 1
#                 label_color.append(n)
#     return label_color

def volcano_plot(deseq_table, lfc_cutoff=0.5, pvalue_cutoff=1e-5, figsize=(4,4), add_text='auto', n_gene_name=15, plot_gene_names = [],
                 fig_title='', save_path='volcano_plot.pdf'):
    '''
    add_text: 'auto' or 'define' or 'none'
    '''
    result = pd.DataFrame()
    result['x'] = deseq_table['log2FoldChange']
    min_value_excluding_0 = np.min(deseq_table[deseq_table['padj'] != 0]['padj'])
    deseq_table.loc[deseq_table['padj'] == 0, 'padj'] = min_value_excluding_0
    result['y'] = -np.log10(deseq_table['padj'])
    # max_value_exclueding_inf = np.max(result[result['y'] != np.inf]['y'])
    # result.loc[result['y'] == np.inf, 'y'] = max_value_exclueding_inf + 50

    # set cutoff
    x_threshold = lfc_cutoff
    y_threshold = -np.log10(pvalue_cutoff)
    
    # group to up, normal, down
    result['group'] = 'grey'
    result.loc[(result.x > x_threshold)&(result.y > y_threshold),'group'] = 'tab:red' # up regulated
    result.loc[(result.x < -x_threshold)&(result.y > y_threshold),'group'] = 'tab:blue' # down regulated
    # print(result.sort_values('y', ascending=False).head(10))
    # set which gene use to plot gene name
    if add_text == 'auto':
        result['xy'] = result['x'] * result['y']
        up_gene_names = result[result['group'] == 'tab:red'].sort_values('xy', ascending=False).head(n_gene_name//2).index
        down_gene_names = result[result['group'] == 'tab:blue'].sort_values('xy', ascending=True).head(n_gene_name//2).index
        plot_gene_names = list(up_gene_names) + list(down_gene_names)
        # plot_gene_names = result[result['group'] != 'grey'].sort_values('xy', ascending=False).head(n_gene_name).index
    # print(result[result['group'] != 'grey'].sort_values('xy', ascending=True).head(n_gene_name//2))
    # print(up_gene_names)
    # print(down_gene_names)
    # plot
    fig, ax = plt.subplots(figsize=figsize, dpi=600) 
    ax.set(title=fig_title)
    ax.scatter(result['x'], result['y'], s=2, alpha=.5, c=result['group'])
    ax.set_ylabel('-log10(P adj.)',fontweight='bold', fontsize=10)
    ax.set_xlabel('log2(Fold change)',fontweight='bold', fontsize=10)
    ax.spines['right'].set_visible(False) 
    ax.spines['top'].set_visible(False) 

    plt.axvline(-x_threshold, color='black',linestyle='dashed', linewidth=.5) 
    plt.axvline(x_threshold, color='black',linestyle='dashed', linewidth=.5) 
    plt.axhline(y_threshold, color='black',linestyle='dashed', linewidth=.5) 
    if add_text == 'auto' or add_text == 'define':
        get_name_color = lambda txt: '#FD693E' if txt in up_gene_names else '#1873B8' if txt in down_gene_names else '#000000'
        texts = [ax.annotate(txt, (result.loc[txt, 'x'], result.loc[txt, 'y']),
                             color=get_name_color(txt), fontsize=6, fontweight='bold') for txt in plot_gene_names]
        # print(texts)
        adjustText.adjust_text(texts, 
                               arrowprops=dict(arrowstyle='-', color='black', lw=0.2),
                               expand=(1.2, 1.2),
                               force_text=(0.2, 0.5),
                               force_explode=(0.2, 2),
                               avoid_self=True,
                               prevent_crossings=True,
                               ensure_inside_axes=True,
                               expand_axes=False,
                               max_move=None)
    fig.tight_layout()
    fig.savefig(save_path, bbox_inches = 'tight')


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
                 figsize=(3,3), 
                 add_text=add_text, 
                 n_gene_name=n_gene_name, 
                 plot_gene_names = plot_gene_names,
                 fig_title=fig_title, 
                 save_path=save_path)

if __name__ == '__main__':
    main()