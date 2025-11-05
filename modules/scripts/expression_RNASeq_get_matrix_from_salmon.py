#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   expression_RNASeq_get_matrix_from_salmon.py
@Time    :   2025/11/05 11:02:33
@Author  :   Xin Dong
@Contact :   xindongnov@gmail.com
@License :   (C)Copyright 2020-2025, XinDong
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
    parser.add_argument('-i', '--input', help='input files', required=True)
    parser.add_argument('-g', '--gene_table', help='the gene table path for mapping gene name', default=None)
    parser.add_argument('-c', '--count', help='the count path to save', default="./count.txt")
    parser.add_argument('-t', '--tpm', help='the TPM path to save', default="./tpm.txt")
    args = parser.parse_args()

    files = args.input.split(',')
    # input_type = args.format
    count_path = args.count
    tpm_path = args.tpm
    

    df_list = [pd.read_csv(i, sep='\t', index_col=0, comment='#') for i in files]
    names = [i.split('/')[2] for i in files]

    count_ls = []
    tpm_ls = []

    # if input_type == 'salmon':
    for i, l in enumerate(df_list):
        count_tmp = l.rename(columns = {'NumReads' : names[i]})
        count_ls.append(count_tmp[names[i]])
        tpm_tmp = l.rename(columns = {'TPM' : names[i]})
        tpm_ls.append(tpm_tmp[names[i]])
    # elif input_type == 'STAR':
        # pass
        # for i, l in enumerate(df_list):
        #     count_tmp = l.rename(columns = {'NumReads' : names[i]})
        #     count_ls.append(count_tmp[names[i]])
        #     tpm_tmp = l.rename(columns = {'TPM' : names[i]})
        #     tpm_ls.append(tpm_tmp[names[i]])
    # else:
        # sys.exit(0)

    count = pd.concat(count_ls, axis=1)
    tpm = pd.concat(tpm_ls, axis=1)

    count.to_csv(count_path, sep='\t', index=True, header=True)
    tpm.to_csv(tpm_path, sep='\t', index=True, header=True)


if __name__ == "__main__":
    main()