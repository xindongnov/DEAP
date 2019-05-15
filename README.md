# DEAP

## Introducton

Differential Expression Analysis Pipeline(DEAP) is a pipeline which can handle expression data that RNA-seq or MicroArray experiment generated. The pipeline is built using [snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) which allows for ease of use, optimal speed, and a highly modular code that can be further added onto and customized by experienced users. 

## Usage
```{bash}
conda env create -n DEAP -f environment.yaml
conda activate DEAP
```
Then make your folder structure like this:  
 
data  *(store your data)  
DEAP *(this one is the git repo)*  
metasheet.csv  
ref.yaml  
config.yaml  

At last type:
```
snakemake -s DEAP/DEAP.snakefile
```

## Appendix
