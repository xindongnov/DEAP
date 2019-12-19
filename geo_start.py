#!/usr/bin/env python

# ================================
# @auther: Xin Dong
# @email: xindong9511@gmail.com
# @date: Dec 2019
# ================================

import os
import sys
import subprocess
import pandas as pd
import yaml
import argparse


def main():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser()
    parser.add_argument('-s', '--sampleconfig', help='The path of new config file with samples.', default="sample_config.yaml")
    parser.add_argument('-e', '--command', help='Other commands you want to submit to snakemake, should be quoted. This would overwrite "-j"/"--jobs". Eg. "-npr"', default="")
    parser.add_argument('-j', '--jobs', help='The cores you want to provide to snakemake', default="8")
    args = parser.parse_args()

    metadata = pd.read_csv("metasheet.csv", index_col=0, sep=',', comment='#', skipinitialspace=True)
    try:
        os.mkdir("data")
    except FileExistsError:
        sys.stdout.write("Folder data has been created!\n")
    download_flag = True
    # laytypedic = {}
    id_run_list = list(zip(metadata["sample"],metadata.index))
    failed_list = []
    for i,j in id_run_list:
        try:
            sys.stdout.write("Downloading sample %s \n" % i)
            ret = subprocess.run("python modules/scripts/geo_rawdown.py -i {id} -o {path} -g".format(id=i,path="data/%s" % j), 
                                shell=True, check=True)
        except subprocess.CalledProcessError:
            sys.stderr.write("Sample %s failed to download!\n" % i)
            failed_list.append(i)
            download_flag = False

    if download_flag == True:
        sys.stdout.write("Generating new config.yaml file. Store to %s." % args.sampleconfig)
        f = open("config.yaml")
        config = yaml.safe_load(f)
        f.close()
        config = {}
        config["samples"] = {}
        for i,j in id_run_list:
            folder = os.listdir("data/%s" % j)
            config["samples"][i] = [os.path.join("data/%s" % j,file_name) for file_name in folder if i in file_name]
        with open(args.sampleconfig,"w") as sample_config:
            yaml.safe_dump(config, sample_config)
        snkmk_cmd = "snakemake -s DEAP/DEAP.snakefile --configfile %s -j %s %s" % (args.sampleconfig,args.jobs,args.command)
        os.system(snkmk_cmd)
    else:
        sys.stderr.write("Following samples failed in download. Please manually check these samples:\n")
        sys.stderr.write("\t".join(failed_list))

if __name__ == "__main__":
    main()


    