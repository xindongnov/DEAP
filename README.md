# DEAP (Alpha)

## Introducton

Differential Expression Analysis Pipeline(DEAP) is a pipeline that can handle expression data that RNA-seq or MicroArray experiments generated. The pipeline is built using [snakemake](https://snakemake.readthedocs.io/), which allows for ease of use, optimal speed, and highly modular code that can be further added onto and customized by experienced users.

## Basic Usage

We suppose you start your project with `PROJECT` folder:

```bash
mkdir PROJECT
cd PROJECT
git clone git@github.com:XinDong9511/DEAP.git
conda env create -n DEAP -f DEAP/environment.yaml # this might took long time based on your network environment
conda activate DEAP
```

The data can be store anywhere with given path in config.yaml file, but the best practice is storing all your data in `data` folder.

```bash
mkdir data
```

You can download reference files [/mnt/Storage/home/dongxin/Files/DEAP_ref_files](blank). All compressed files should be extracted to `ref_files` folder.

```bash
mkdir ref_files
cd ref_files
wget https://XXXXXXX/hg38.tar.gz
wget https://XXXXXXX/mm10.tar.gz
wget https://XXXXXXX/GPL.tar.gz
tar xvzf hg38.tar.gz
tar xvzf mm10.tar.gz
tar xvzf GPL.tar.gz
cd ..
```

If you have already downloaded reference files on your server or local machine, we recommend you link the reference files by soft link:

```bash
ln -s /mnt/Storage/home/dongxin/Files/DEAP_ref_files ref_files
```

Copy needed config files to current folder.

```bash
cp DEAP/metasheet.csv .
cp DEAP/ref.yaml .
cp DEAP/config.yaml .
```

Then your folder structure should like this:  

> PROJECT/  
├── data  *(store your data)  
├── ref_files  
├── DEAP *(this one is the git repo)*  
├── metasheet.csv  
├── ref.yaml  
└── config.yaml  

Please modify your config.yaml and metasheet following the intruction in **Files Format** part. After modified yaml files and metasheet.csv, at last type:

```bash
snakemake -s DEAP/DEAP.snakefile
```

Then the pipeline would start running for analysis. If you only want to check the tasks scheduling and command used in analysis, you can type:

```bash
snakemake -s DEAP/DEAP.snakefile -npr
```

Other snakemake advanced instruction can be found [here](https://snakemake.readthedocs.io/).

## All Data From GEO

[Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/gds/) database stores a huge number curated gene expression DataSets, as well as original Series and Platform records. If you want to process data all from GEO, you can use the `geo_start.py` to help you download data and configure pipeline.

```shell
usage: geo_start.py [-h] [-s SAMPLECONFIG] [-m METASHEET] [-e COMMAND]
                    [-j JOBS]

optional arguments:
  -h, --help            show this help message and exit
  -s SAMPLECONFIG, --sampleconfig SAMPLECONFIG
                        The path of new config file with samples.
  -e COMMAND, --command COMMAND
                        Other command you want to submit to snakemake,
                        should be quoted. This would overwrite "-j"/
                        "--jobs". Eg. "-npr"
  -j JOBS, --jobs JOBS  The cores you want to provide to snakemake
```

Basic configuration is same as mentioned in **Basic Usage** part.

```bash
mkdir PROJECT
cd PROJECT
git clone git@github.com:XinDong9511/DEAP.git
conda env create -n DEAP -f DEAP/environment.yaml # this might took long time based on your network environment
conda activate DEAP

# use this
mkdir ref_files
cd ref_files
wget https://XXXXXXX/hg38.tar.gz
wget https://XXXXXXX/mm10.tar.gz
wget https://XXXXXXX/GPL.tar.gz
tar xvzf hg38.tar.gz
tar xvzf mm10.tar.gz
tar xvzf GPL.tar.gz
cd ..
# or this
ln -s path/to/previous/ref_files/ ref_files

cp DEAP/ref.yaml .
cp DEAP/metasheet.csv .
cp DEAP/config.yaml .
```

In this case, you still need to modify the metasheet table to define samples relationship as instructed from **Files Format** first. Please noticed that all sample ID (2nd columns) should be exact accession with upper case "GSM", eg. `GSM123456`.  
But in config.yaml file, you can skip the `samples` part just leave them as blank.  

After configuration of these two files, you can start download and processing by:

```bash
python DEAP/geo_start.py
```

This script also support passing parameter to snakemake:

```bash
python DEAP/geo_start.py -j 8 # use 8 core to run snakemake
python DEAP/geo_start.py -e "-npr" # dry-run for snakemake
```

You also can specify the output new config files:

```bash
python DEAP/geo_start.py -s new_config.yaml # dry-run for snakemake
```

## Files Format

- config.yaml:  
Here are several keys in this yaml file.  
  - `metasheet`: a path for a table file that defined the relationship among samples.  
  - `ref`: a yaml file that stored the path where reference files are.  
  - `aligner`: the align tool that you want to use for RNA-seq samples, options: `"STAR"` or `"salmon"`.  
  - `assembly`: the assembly you want to use for your data, options: `"mm10"` or `"hg38"`.  
  - `trim`: DEAP use [trim_galore](https://github.com/FelixKrueger/TrimGalore) to cut off the adaper at the distal of reads and do quality control, options: `True` or `False`.  
  - `lisa`: [Lisa](http://lisa.cistrome.org/) is used to determine the transcription factors and chromatin regulators that are directly responsible for the perturbation of a differentially expressed gene set. option:`True` or `False`.  
  - `check_compare`: options: `True` or `False`, when set `True` for this option, DEAP will check comparison relationship and won't process alignment for RNA-seq samples that do not have enough samples(at least 2 samples in one condition). When set `False`, DEAP will continue processing trimming adaptor and alignment but still won't do comparison.
  - `samples`: Can include many values of samples and in each key, you could write one fastq file for single-end smaples or two fastq files for pair-end samples. *You can not add this part if you use `geo_start.py` to run pipeline.*  

- metasheet.csv:  
This is a table that defined samples relationship and define how you would like to do differential expression genes among samples. Except `compare_XXX`, other headers should not be modified.  
  - **The 1st column (run)** is the run name. A run(or batch) of samples should share same name though they could be devided into several rows for different FASTQ files or CEL files.  
  - **The 2nd column (sample)** should always be samples ID. This sample ID **must** exactly match the sample names configured in the config yaml file. *If you use `geo_start.py` to run pipeline, please make sure all sample ID is exact accession with upper case "GSM", eg. `GSM123456`.*  
  - **The 3rd column (experment_type)** defined the experiment type of data, please choose one of `{RS, MA_A, MA_O, MA_I}`. `RS` for RNA-seq, `MA_{}` for microarry data generated by **A**ffymatrix, **O**ligo or **I**llumina *(Not support yet)*.  
  - **The 4th colunm (platform)**: You are supposed to write the platform of each sample. We only use the first one platform you defined in same run but still recommand you fulfill the whole table. For RNA-seq, the platform will be used for reads group. For microarray data, the platform should use GPL file you downloaded from GEO, eg. `GPL570`. Please make sure the GPL you write in table, has the matched file which named by `GPL{XXXX}.txt` in `ref_files/GPL/`. The `{XXXX}` is the number of this platform.
  - **The 5th colunm (treatment)**: The annotation of each sample. You should write some string in this column. For example, `"Control"`, `"RNAi"`, `"Treat"`, `"ESR1i"` etc.
  - **The 6th column (compare_XXX)**: The samples that you want to perform a Differential Expression on using limma and DEseq. The “control” should be marked with a `1`, and the “treatment” should be marked with a `2`. You are allowed to have as many `compare_XXXX` columns as you want at the end.  
  
- ref.yaml:  
Here you are allowed to define your own reference files. Basically includes GPL files, hg38 and mm10 index files .  

## Appendix

### LISA installation

LISA isn't a necessary component of DEAP, you needn't install is if you do not plan to run it -- set `lisa` in config as `False`. Otherwise, please install LISA by following instruction.  

Since LISA also uses snakemake to build the workflow, to avoid the version conflict, we recommend you built a new environment for LISA.  
[LISA GIT](https://github.com/qinqian/lisa)

1. Install LISA by conda:

    ```bash
    conda create -n lisa python=3.6
    conda activate lisa
    conda install -c qinqian lisa
    ```

2. Download reference files of LISA:

    ```bash
    wget --user=lisa --password='xxx'  http://lisa.cistrome.org/cistromedb_data/lisa_v1.0_hg38.tar.gz
    wget --user=lisa --password='xxx'  http://lisa.cistrome.org/cistromedb_data/lisa_v1.1_mm10.tar.gz
    ```

3. Extract files and update LISA configuration: *The files can not move after `lisa_update_conf`, please extract them to a proper position.*

    ```bash
    tar xvfz lisa_v1.0_hg38.tar.gz
    lisa_update_conf --folder hg38/ --species hg38
    # or
    tar xvfz lisa_v1.0_mm10.tar.gz
    lisa_update_conf --folder mm10/ --species mm10
    ```

4. Exit LISA environment and activate DEAP:

    ```bash
    conda deactivate
    conda activate DEAP
    ```
