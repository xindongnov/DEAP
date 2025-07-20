#!/usr/bin/env python

# =================
# @File    :   DEAP.snakefile
# @Time    :   2023/07/14 21:27:34
# @Author  :   Xin Dong
# @Contact :   xindong9511@gmail.com
# @License :   (C)Copyright 2020-2023, XinDong
# =================


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
    _invalid_map = {'\r':'\n', '$':''}
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
#             comps_info[comp[8:]]['control'] = meta_info[meta_info[comp] == 0].index
#             comps_info[comp[8:]]['treat'] = meta_info[meta_info[comp] == 1].index
#     return comps_info

def updateMeta(config):
    # _sanity_checks(config)
    metadata = pd.read_csv(config['metasheet'], index_col=0, sep=',', comment='#', skipinitialspace=True)
    config["runs"] = {}
    for run in set(metadata.index.values):
        if isinstance(metadata.loc[run,"experiment_type"],str):
            sys.stderr.write("ERROR: %s does NOT match any mates." % run)
            sys.exit(1)
        elif list(metadata.loc[run,"experiment_type"])[0] in ["RS","MA_A","MA_O"]:
            config["runs"][run] = {'type': list(metadata.loc[run,"experiment_type"])[0], 'platform': list(metadata.loc[run,"platform"])[0]}
            comp_list = ['sample', 'condition', 'batch']
            comp_list.extend(metadata.columns.values[5:]) # 5 means after 'batch'
            config["runs"][run]['samples'] = {}
            config["runs"][run]['compare'] = {}
            design = metadata.loc[run,comp_list]
            for s in list(metadata.loc[run,"sample"]):
                config["runs"][run]['samples'][s] = config['samples'][s]
            for c in design.loc[run,comp_list[3:]]: # 0 is sample, 1 is condition, 2 is batch
                if len(pd.unique(design.loc[:,c])) != 1:
                    config["runs"][run]['compare'][c] = {'control': {'name':list(design.loc[design.loc[:,c] == 0,'condition'])[0], 
                                                                        'sample': list(design.loc[design.loc[:,c] == 0,'sample'])},
                                                            'treat': {'name':list(design.loc[design.loc[:,c] == 1,'condition'])[0], 
                                                                        'sample': list(design.loc[design.loc[:,c] == 1,'sample'])}}
            config["runs"][run]['raw_design'] = metadata.loc[run,comp_list]
        else:
            sys.stdout.write("WARNING: %s does NOT match any Experiment type." % run)
    return config


# def add_lisa_config(config):
#     conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
#     config["conda_root"] = conda_root
#     conda_path = os.path.join(conda_root, 'pkgs')
#     if not "lisa_path" in config or not config["lisa_path"]:
#         config["lisa_path"] = os.path.join(conda_root, 'envs', 'lisa', 'bin', 'lisa')
#         config["lisa_env"] = os.path.join(conda_root, 'envs', 'lisa', 'bin')


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
# add_lisa_config(config)

#NOW load ref.yaml - SIDE-EFFECT: loadRef CHANGES config
loadRef(config)
#-----------------------------------------

global RES_PATH
RES_PATH = config['res_path']

def all_targets(wildcards):
    print(config)
    ls = []
    #IMPORT all of the module targets
    if config['trim'] == True:
        ls.extend(trim_targets(wildcards))
    if config['aligner'] == 'STAR':
        ls.extend(align_STAR_targets(wildcards))
        # ls.extend(rsem_quantification_targets(wildcards))
    else:
        ls.extend(align_salmon_targets(wildcards))
    ls.extend(experssion_targets(wildcards))
    if config['lisa'] == True:
        ls.extend(lisa_targets(wildcards))
    ls.extend(geneontology_targets(wildcards))
    return ls

rule all:
    input: all_targets

include: "./modules/trim.snakefile"
if config['aligner'] == 'STAR':
    include: "./modules/align_STAR.snakefile"  # rules specific to STAR
    # include: "./modules/rsem.snakefile"        # rules specific to STAR
else:
    include: "./modules/align_salmon.snakefile"   # rules specific to salmon
include: "./modules/expression.snakefile"
if config['lisa'] == True:
    include: "./modules/lisa.snakefile"
include: "./modules/geneontology.snakefile"  # rules specific to gene ontology

