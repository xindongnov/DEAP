# align salmon 


def align_salmon_targets(wildcards):
    ls = []
    for run in config["RS_runs"]:
        for sample in run['sample']:
        	ls.append('analysis/%s/%s/')
    return ls

def getFastq(wildcards):
    return config[]


rule align_salmon:
    input:
        ""
    output:
        ""
    params:


