#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   expression_remove_batch.py
@Time    :   2024/08/19 19:31:30
@Author  :   Xin Dong
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2024, XinDong
'''

import sys
from inmoose.pycombat import pycombat_norm, pycombat_seq
import argparse

def main():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser()
    parser.add_argument('-d', '--design', help='The design file', required=True)
    parser.add_argument('-c', '--count', help='the raw count table path', default="count.txt")
    parser.add_argument('-cmp', '--compare', help='the compare type', default='treat_vs_ref')
    parser.add_argument('-l', '--lfc', help='the lfc cutoff', default=1.5)
    parser.add_argument('-p', '--pvalue', help='the pvalue cutoff', default=1e-5)
    parser.add_argument('-o', '--output', help='the output path', default='DEseq.txt')
    args = parser.parse_args()
    
    # remove batch effect
    for i in metadata['batch'].unique():
        if metadata[metadata['batch'] == i].__len__() == 1:
            do_batch=False
    if do_batch:
        counts_df = pycombat_seq(counts_df.T,metadata['batch']).T