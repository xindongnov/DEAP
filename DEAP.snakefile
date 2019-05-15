#!/usr/bin/env python

import os
import sys
import subprocess
import numpy as np
import pandas as pd
import yaml
from string import Template
from collections import defaultdict
# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

def _sanity_checks(config):
    #metasheet pre-parser: converts dos2unix, catches invalid chars
    _invalid_map = {'\r':'\n', '(':'.', ')':'.', '$':''}
    _meta_f = open(config['metasheet'])
    _meta = _meta_f.read()
    _meta_f.close()

    _tmp = _meta.replace('\r\n','\n')
    #check other invalids
    for k in _invalid_map.keys():
        if k in _tmp:
            _tmp = _tmp.replace(k, _invalid_map[k])

    #did the contents change?--rewrite the metafile
    if _meta != _tmp:
        #print('converting')
        _meta_f = open(config['metasheet'], 'w')
        _meta_f.write(_tmp)
        _meta_f.close()

# def _get_comp_info(meta_info):
#     comps_info = defaultdict(dict)
#     for comp in meta_info.columns:
#         print(meta_info)
#         if comp[:8] == 'compare_':
#             comps_info[comp[8:]]['control'] = meta_info[meta_info[comp] == 1].index
#             comps_info[comp[8:]]['treat'] = meta_info[meta_info[comp] == 2].index
#     return comps_info

def updateMeta(config):
    _sanity_checks(config)
    metadata = pd.read_csv(config['metasheet'], index_col=0, sep=',', comment='#', skipinitialspace=True)
    config["RS_runs"] = {}
    config["MA_runs"] = {}
    for run in set(metadata.index.values):
        if isinstance(metadata.loc[run,"experment_type"],str):
            if metadata.loc[run,"experment_type"].startswith("MA_"):
                config["MA_runs"][run] = {'type': metadata.loc[run,"experment_type"],
                                          'samples': config['samples'][metadata.loc[run,"sample"]],
                                          'matrix': metadata.loc[run,"condition"]}
            elif metadata.loc[run,"experment_type"] == "RS":
                sys.stdout.write("ERROR: %s does NOT match any mates." % run)
                sys.exit(1)
            else:
                sys.stdout.write("WARNING: %s does NOT match any Experiment type." % run)
        else:
            if list(metadata.loc[run,"experment_type"])[0] == "RS":
                config["RS_runs"][run] = {'type': 'RS'}
                # wrote sample information into config
                # config["RS_runs"][run]['samples'] = {}
                # for s in list(metadata.loc[run,"sample"]):
                #     config["RS_runs"][run]['samples'][s] = config['samples'][s]
                # wrote design list into config
                comp_list = ['sample','treatment']
                comp_list.extend([i for i in metadata.columns.values if i.startswith('compare_')])
                # handle none replicate
                design = metadata.loc[run,comp_list]
                # will not run unlist sample
                config["RS_runs"][run]['samples'] = {}
                for c in design.loc[run,comp_list[2:]]:
                    # print(design.loc[:,c])
                    sample_list = []
                    control = list(design.loc[:,c]).count(1)
                    treat = list(design.loc[:,c]).count(2)
                    NA = len(design.loc[:,c]) - treat - control
                    if control == 1 or treat == 1 or NA == len(design.loc[:,c]):
                        sys.stdout.write("WARNING: run: %s, compare: %s do not have enough data!\n" % (run,c))
                        sys.stdout.write("The samples in this comparison will not implement alignment and differential expression.\n")
                    else:
                        sample_list.extend(design.loc[design.loc[:,c] == 1,'sample'])
                        sample_list.extend(design.loc[design.loc[:,c] == 2,'sample'])
                        # print(sample_list)
                        for s in sample_list:
                            if s not in config["RS_runs"][run]['samples']:
                                config["RS_runs"][run]['samples'][s] = config['samples'][s]
                        # print(config["RS_runs"][run]['samples']
                        # print(design.loc[design.loc[:,c].notna(),'sample'])
                config["RS_runs"][run]['compare'] = design
            elif metadata.loc[run,"experment_type"].startswith("MA_"):
                sys.stdout.write("ERROR: %s has more than one microarray folder." % run)
                sys.exit(2)
            else:
                sys.stdout.write("WARNING: %s does NOT match any Experiment type." % run)
    # print(config)
    # config["comparisons"] = [c[8:] for c in metadata.columns if c.startswith("compare_")]
    # config["comps"] = _get_comp_info(metadata)
    # config["metacols"] = [c for c in metadata.columns if c.lower()[:4] != 'compare']
    # config["file_info"] = { sampleName : config["samples"][sampleName] for sampleName in metadata.index }
    return config

def loadRef(config):
    """Adds the static reference paths found in config['ref']
    NOTE: if the elm is already defined, then we DO NOT clobber the value
    """
    f = open(config['ref'])
    ref_info = yaml.safe_load(f)
    f.close()
    #print(ref_info[config['assembly']])
    for (k,v) in ref_info[config['assembly']].items():
        #NO CLOBBERING what is user-defined!
        if k not in config:
            config[k] = v

#---------  CONFIG set up  ---------------
configfile: "config.yaml"   # This makes snakemake load up yaml into config 
config = updateMeta(config)
# addPy2Paths_Config(config)

#NOW load ref.yaml - SIDE-EFFECT: loadRef CHANGES config
loadRef(config)
#-----------------------------------------

def all_targets(wildcards):
    # print(config)
    ls = []
    #IMPORT all of the module targets
    if config['trim'] == True:
        ls.extend(trim_targets(wildcards))
    ls.extend(align_salmon_targets(wildcards))
    ls.extend(DE_RNAseq_targets(wildcards))

    return ls   

rule all:
    input: all_targets

include: "./modules/trim.snakefile"

if config['aligner'] == 'STAR':
    include: "./modules/align_STAR.snakefile"     # rules specific to STAR
else:
    include: "./modules/align_salmon.snakefile"   # rules specific to salmon

include: "./modules/DE_RNAseq.snakefile"



