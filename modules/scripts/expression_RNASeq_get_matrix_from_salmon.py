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
    parser.add_argument('-i', '--input', help='input sf files', required=True)
    parser.add_argument('-o', '--output', help='the path to save', default="./output.txt")
    parser.add_argument('-n', '--name', help='input names', required=True)
    parser.add_argument('-t', '--type', help='type of data', choices = ['Rawcount', 'TPM'], required=True)
    args = parser.parse_args()

    files = args.input.split(',')
    result_path = args.output
    names = args.name.split(',')
    dtype = args.type

    df_list = [pd.read_csv(i, sep='\t', index_col=0) for i in files]
    
    res = []
    if dtype == 'Rawcount':
        for i, l in enumerate(df_list):
            tmp = l.rename(columns = {'NumReads' : names[i]})
            res.append(tmp[names[i]])
    elif dtype == 'TPM':
        for i, l in enumerate(df_list):
            tmp = l.rename(columns = {'TPM' : names[i]})
            res.append(tmp[names[i]])

    output = pd.concat(res, axis=1)

    output.to_csv(result_path, sep='\t', index=True, header=True)


if __name__ == "__main__":
    main()