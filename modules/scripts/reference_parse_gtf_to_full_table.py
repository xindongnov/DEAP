#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   reference_prepare.py
@Time    :   2023/07/14 16:48:41
@Author  :   Xin Dong
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2023, XinDong
'''

# this script is used to parse the gtf file to a full table, also output is gzipped
# usage: python reference_prepare.py -g <gtf file> -o <output file>


import argparse

def read_gtf(gtf_path):
    import pandas as pd
    gtf = pd.read_csv(gtf_path, sep='\t', comment='#', header=None, low_memory=False)
    print('Parsing information ...')
    anno_dict = {}
    for n, i in enumerate(gtf[8]):
        l = [item.split(' "') for item in i.split('; ')]
        anno_dict[n] = {}
        for k, j in l:
            anno_dict[n][k] = j.strip(';').strip('"')
    print('Post-processing ...')
    coor = gtf.loc[:,0:7].copy()
    coor.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame']
    ddf = pd.DataFrame(anno_dict).T
    res = pd.concat([coor, ddf], axis=1)
    res = res.fillna('NA')
    print('Done!')
    return res

if __name__ == '__main__':      
    args = argparse.ArgumentParser()
    args.add_argument('-g', '--gtf', help='raw gtf files', type=str, required=True)
    args.add_argument('-o', '--output', help='output txt files', type=str, required=True)
    args = args.parse_args()

    gtf = args.gtf
    output = args.output
    if not output.endswith('.gz'):
        output = output + '.gz'

    read_gtf(gtf).to_csv(output, sep='\t', index=False, compression='gzip')


    