#!/usr/bin/env python

import os
import sys
import subprocess
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -l params.lisa_path -s params.species -p params.prefix -t threads -i input"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-l", "--lisa_path", help="lisa_path")
    optparser.add_option("-s", "--species", help="species")
    optparser.add_option("-p", "--prefix", help="prefix")
    optparser.add_option("-t", "--threads", help="threads")
    optparser.add_option("-i", "--input", help="input")
    optparser.add_option("-c", "--conda", help="conda_root")
    (options, args) = optparser.parse_args(sys.argv)

    cmd = ". %s/etc/profile.d/conda.sh && conda activate lisa; " % options.conda
    cmd += "%s model --method=\"all\" --web=False --new_rp_h5=None --new_count_h5=None --species %s --epigenome \"['DNase', 'H3K27ac']\" --cluster=False --covariates=False --random=True --prefix %s --background=dynamic_auto_tad --stat_background_number=1000 --threads %s %s " % (options.lisa_path,options.species,options.prefix,options.threads,options.input)
    cmd += "&& conda deactivate" #% options.conda
    print(cmd)
    subprocess.run(cmd,shell=True,executable='/bin/bash')

if __name__ == '__main__':
    main()

