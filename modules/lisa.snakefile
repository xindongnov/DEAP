
# ================================
# @auther: Xin Dong
# @email: xindong9511@gmail.com
# @date: Sep 2019
# ================================


def lisa_targets(wildcards):
    print(snakemake.snakemake)
    ls = []
    for run in config["runs"]:
        for comp in config["runs"][run]['compare']:
            ctrl = config["runs"][run]['compare'][comp]['control']['name']
            treat = config["runs"][run]['compare'][comp]['treat']['name']
            if len(config["runs"][run]['compare'][comp]['control']['sample']) > 1 and len(config["runs"][run]['compare'][comp]['treat']['sample']) > 1:
                for geneType in ["upRegGenes","downRegGenes"]:
                    ls.append("analysis/%s/lisa/%s_%s/%s_%s.%s.txt" % (run,comp,geneType,run,comp,geneType))
                    # n = subprocess.check_output("wc -l %s" % checkpoints.lisa_copy_file.get(run=run,compare=comp,geneType=geneType).output[0],shell=True)
                    # n = int(n.decode("utf-8").strip().split()[0])
                    # if n != 0:
                    ls.append("analysis/%s/lisa/%s_%s/%s_%s.%s.txt.1000.lisa_direct.csv" % (run,comp,geneType,run,comp,geneType))
    # for run in config["RS_runs"]:
    #     for comp in config["RS_runs"][run]['compare']:
    #         ctrl = config["RS_runs"][run]['compare'][comp]['control']['name']
    #         treat = config["RS_runs"][run]['compare'][comp]['treat']['name']
    #         if len(config["RS_runs"][run]['compare'][comp]['control']['sample']) > 1 and len(config["RS_runs"][run]['compare'][comp]['treat']['sample']) > 1:
    #             for geneType in ["upRegGenes","downRegGenes"]:
    #                 ls.append("analysis/%s/lisa/%s_%s_%s/%s_%s.%s.txt" % (run,treat,ctrl,geneType,treat,ctrl,geneType))
    #                 n = subprocess.check_output("wc -l %s" % checkpoints.lisa_copy_file.get(run=run,treatment="%s_%s" % (treat,ctrl),geneType=geneType).output[0],shell=True)
    #                 n = int(n.decode("utf-8").strip().split()[0])
    #                 if n != 0:
    #                     ls.append("analysis/%s/lisa/%s_%s_%s/%s_%s.%s.txt.1000.lisa_direct.csv" % (run,treat,ctrl,geneType,treat,ctrl,geneType))
                # ls.append("analysis/%s/lisa/%s_%s_upRegGenes/%s_%s.upRegGenes.txt" % (run,treat,ctrl,treat,ctrl))
                # ls.append("analysis/%s/lisa/%s_%s_downRegGenes/%s_%s.downRegGenes.txt" % (run,treat,ctrl,treat,ctrl))
                # n = subprocess.check_output("wc -l %s" % checkpoints.lisa_copy_file.get(run=run,treatment="%s_%s" % (treat,ctrl),geneType="downRegGenes").output[0],shell=True).decode("utf-8").strip().split()[0]
                # n = int(n)
                # if n != 0:
                #     ls.append("analysis/%s/lisa/%s_%s_downRegGenes/%s_%s.downRegGenes.txt.1000.lisa_direct.csv" % (run,treat,ctrl,treat,ctrl))
                # ls.append("analysis/%s/lisa/%s_%s_upRegGenes/%s_%s.upRegGenes.txt.1000.lisa_direct.csv" % (run,treat,ctrl,treat,ctrl))
    return ls


rule lisa_copy_file:
    input:
        "analysis/{run}/expression/{compare}/{run}_{compare}.{geneType}.txt"
    output:
        "analysis/{run}/lisa/{compare}_{geneType}/{run}_{compare}.{geneType}.txt"
    shell:
        "head -n 500 {input} > {output}"

# rule lisa_copy_down_file:
#     input:
#         "analysis/{run}/expression/{run}_results/{treatment}.{geneType}.txt"
#     output:
#         "analysis/{run}/lisa/{treatment}_{geneType}/{treatment}.{geneType}.txt"
#     shell:
#         "cp {input} {output}"

rule lisa_run:
    input:
        "analysis/{run}/lisa/{compare}_{geneType}/{run}_{compare}.{geneType}.txt"
    output:
        "analysis/{run}/lisa/{compare}_{geneType}/{run}_{compare}.{geneType}.txt.1000.lisa_direct.csv",
        "analysis/{run}/lisa/{compare}_{geneType}/{run}_{compare}.{geneType}.txt.profile",
        "analysis/{run}/lisa/{compare}_{geneType}/{run}_{compare}.{geneType}.txt.Snakefile.model",
        "analysis/{run}/lisa/{compare}_{geneType}/{run}_{compare}.{geneType}.txt.yml"
    params:
        prefix=lambda wildcards: "%s_%s.%s.txt" % (wildcards.run,wildcards.compare,wildcards.geneType),
        species="hg38",
        lisa_path=config["lisa_path"],
    threads: 8
    run:
        cmd = "cd analysis/%s/lisa/%s_%s/; " % (wildcards.run,wildcards.compare,wildcards.geneType)
        cmd += ". %s/etc/profile.d/conda.sh && conda activate lisa; " % config["conda_root"]
        cmd += "%s model --method=\"all\" --web=False --new_rp_h5=None --new_count_h5=None --species %s --epigenome \"['DNase', 'H3K27ac']\" --cluster=False --covariates=False --random=True --prefix %s --background=dynamic_auto_tad --stat_background_number=1000 --threads %s %s; " % (params.lisa_path,params.species,params.prefix,threads,params.prefix)
        cmd += "conda deactivate"
        print(cmd)
        subprocess.run(cmd,shell=True,executable='/bin/bash')
        # # print(cmd2)
        # # subprocess.call(cmd2.split())
        # # print(cmd3)
        # # subprocess.call(cmd3.split())
        """source activate lisa; """
        """{params.lisa_path} model --method="all" --web=False --new_rp_h5=None --new_count_h5=None --species {params.species} --epigenome "['DNase', 'H3K27ac']" """
        """--cluster=False --covariates=False --random=True --prefix {params.prefix} --background=dynamic_auto_tad """
        """--stat_background_number=1000 --threads {threads} {input}; """
        """conda deactivate; """






