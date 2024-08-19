# DEAP (Alpha)

## Introducton

Differential Expression Analysis Pipeline(DEAP) is a pipeline that can process expression data that RNA-seq or MicroArray experiments generated. The pipeline is built using [snakemake](https://snakemake.readthedocs.io/), which allows for ease of use, optimal speed, and highly modular code that can be further added onto and customized by experienced users.

## Basic Usage

We suppose you start your project with `PROJECT` folder:

```bash
mkdir PROJECT
cd PROJECT
git clone git@github.com:xindong95/DEAP.git
conda env create -f DEAP/environment.yaml # this might take long time based on your network environment, use mamba is recommended
conda activate DEAP
```

The data can be stored anywhere with the given path in the config.yaml file, but the best practice is putting all your data in the `data` folder.

```bash
mkdir data
```

You can download reference files [here](/fs/home/dongxin/Files/DEAP_ref_files). All compressed files should be extracted to the `ref_files` folder.

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

If you have already downloaded reference files on your server or local machine, we recommend you link the reference files by the soft link:

```bash
ln -s /fs/home/dongxin/Files/DEAP_ref_files ref_files
```

Copy needed config files to the current folder.

```bash
cp DEAP/metasheet.csv .
cp DEAP/ref.yaml .
cp DEAP/config.yaml .
```

Then your folder structure would be like this:  

> PROJECT/  
├── data  *(store your data)*  
├── ref_files  
├── DEAP *(this one is the git repo)*  
├── metasheet.csv  
├── ref.yaml  
└── config.yaml  

Please modify your config.yaml and metasheet following the instructions in **Files Format** part. After modifing config.yaml and metasheet.csv, type:

```bash
snakemake -s DEAP/DEAP.snakefile
```

Then the pipeline will start running for analysis. If you only want to check the tasks scheduling and command used in the analysis, you can type:

```bash
snakemake -s DEAP/DEAP.snakefile -np
```

Other snakemake advanced instructions can be found [here](https://snakemake.readthedocs.io/).

## All Data Comes From GEO

[Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/gds/) database stores a huge number curated gene expression DataSets, as well as original Series and Platform records. If you want to process data all from GEO, you can use the `geo_start.py` to help you download data and configure the pipeline.

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

The basic configuration is the same as mentioned in **Basic Usage** part.

```bash
mkdir PROJECT
cd PROJECT
git clone git@github.com:xindong95/DEAP.git
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

If you have already downloaded repository and reference files, soft link would be useful for saving space:

```bash
mkdir PROJECT
cd PROJECT
ln -s /path/to/DEAP DEAP
ln -s /path/to/ref_files ref_files
```

In this case, you still need to modify the metasheet table to define samples relationship as instructed from **Files Format** first. Please noticed that all sample ID (2nd columns) should be exact accession with upper case "GSM", eg. `GSM123456`.  
But in config.yaml file, you can skip the `samples` part just leave them as blank.  

After the configuration of these two files, you can start download and processing by:

```bash
python DEAP/geo_start.py
```

This script also supports passing parameter to snakemake:

```bash
python DEAP/geo_start.py -j 8 # use 8 core to run snakemake
python DEAP/geo_start.py -e "-np" # dry-run for snakemake
```

You can specify the output path of new config files as well:

```bash
python DEAP/geo_start.py -s new_config.yaml
```

## Files Format

You can take the config.yaml, metasheet.csv and ref.yaml provided under the DEAP folder as a reference.

- config.yaml:  
Here are several keys in this yaml file.  
  - `metasheet`: a path for a table file that defined the relationship among samples.  
  - `ref`: a yaml file that stored the path where reference files are.  
  - `aligner`: the align tool that you want to use for RNA-seq samples, options: `"salmon"` or `"STAR"` *(Not support yet)*.  
  - `assembly`: the assembly you want to use for your data, options: `"mm10"` or `"hg38"`.  
  - `trim`: DEAP use [trim_galore](https://github.com/FelixKrueger/TrimGalore) to cut off the adapter at the distal of reads and do quality control, options: `True` or `False`.  
  - `batch`: Whether remove batch effect by inmoose.pycombat, options: `True` or `False`.
  - `lisa`: [LISA](http://lisa.cistrome.org/) is used to determine the transcription factors and chromatin regulators that are directly responsible for the perturbation of a differentially expressed gene set. option:`True` or `False`.  
  <!-- - `check_compare`: options: `True` or `False`, when set `True` for this option, DEAP will check comparison relationship and won't process alignment for RNA-seq samples that do not have enough samples(at least two samples in one condition). When set `False`, DEAP will still process trimming adaptor and alignment but still won't make the comparison. -->
  - `samples`: Can include many samples in this key. For each child key, you can write one fastq file for single-end samples or two fastq files for pair-end samples. *You can remove the whole key if you use `geo_start.py` to run the pipeline.*  
  
- metasheet.csv:  
This is a table that define samples' relationship and the way you would like to do differential expression genes among samples. Except `compare_XXX`, other headers should not be modified.  
  - **The 1st column (run)** is the run name. Samples of a run(or batch) should share the same name, though they could be divided into several rows for different FASTQ files or CEL files.  
  - **The 2nd column (sample)** should always be samples ID. This sample ID **must** precisely match the sample names configured in the config yaml file. *If you use `geo_start.py` to run the pipeline, please make sure all sample IDs are exact accession begin with upper-case "GSM", e.g. `GSM123456`.*  
  - **The 3rd column (experment_type)** define the experiment type of data, please choose one of `{RS, MA_A, MA_O, MA_I}`. `RS` for RNA-seq, `MA_{}` for microarray data generated by **A**ffymatrix, **O**ligo or **I**llumina *(Not support yet)*.  
  - **The 4th column (platform)**: You are supposed to write the platform each sample used. We only use the first one platform you defined in the same run but still recommend you fulfill the whole table. For RNA-seq, the platform will be used for the reads group. For microarray data, the platform should use the GPL file you downloaded from GEO, e.g., `GPL570`. Please make sure the GPL you write in the table has the matched file which named by `GPL{XXXX}.txt` in `ref_files/GPL/`. The `{XXXX}` is the number of this platform.
  - **The 5th column (condition)**: The annotation of each sample. It would help if you wrote some string in this column. For example, `"Control"`, `"RNAi"`, `"Treat"`, `"ESR1i"`, etc.
  - **The 6th column (batch)**: The batch information of each sample. Same batch using same number is acceptable.
  - **The 7th column (Any name)**: The samples that you want to perform a Differential Expression to using limma and DEseq. The "control" should be marked with a `0`, and the "treatment" should be marked with a `1`. Recommand using string like `{Treat}_vs_{Control}` for clarity. You are allowed to have as many `{Treat}_vs_{Control}` columns as you want at the end.  
  
- ref.yaml:  
Here you are allowed to define your reference files, includes GPL files, hg38 index, and mm10 index files.  

## Appendix

### LISA installation (updating to lisa2, unfinished)

LISA isn't a necessary component of DEAP, and you needn't install it if you do not plan to run it - set `lisa` in config as `False`. Otherwise, please install LISA based on the following instruction.  

Since LISA also uses snakemake to build the workflow, to avoid the version conflict, we recommend creating a new environment for LISA. Here is the [GIT](https://github.com/qinqian/lisa) of LISA.

1. Install LISA by conda:

    ```bash
    conda create -n lisa python=3.6
    conda activate lisa
    conda install -c qinqian lisa
    ```

2. Download reference files of LISA:

    ```bash
    wget http://lisa.cistrome.org/cistromedb_data/lisa_v1.0_hg38.tar.gz
    wget http://lisa.cistrome.org/cistromedb_data/lisa_v1.1_mm10.tar.gz
    ```

3. Extract files and update LISA configuration: *The files can not move after `lisa_update_conf`, and please extract them to a proper position.*

    ```bash
    tar xvfz lisa_v1.0_hg38.tar.gz
    lisa_update_conf --folder hg38/ --species hg38
    # or
    tar xvfz lisa_v1.0_mm10.tar.gz
    lisa_update_conf --folder mm10/ --species mm10
    ```

4. Exit LISA environment and activate DEAP environment:

    ```bash
    conda deactivate
    conda activate DEAP
    ```

### Reference building

- use script:  
  - salmon_index:  
  `bash DEAP/modules/scripts/reference_salmon_index.sh <reference_output_path> <transcriptome_fasta_path> <threads>`  
  - STAR_index:  
  `bash DEAP/modules/scripts/reference_STAR_index.sh <reference_output_path> <genome_fasta_path> <gtf_path> <threads>`  
  - transformed_gtf:  
  `python DEAP/modules/scripts/reference_prepare.py -g <gtf file> -o <output file>`  

- hg38:
  
  ```yaml
  hg38:
    salmon_index: ./ref_files/hg38/salmon_index/
    STAR_index: ./ref_files/hg38/STAR_index
    gtf: ./ref_files/hg38/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz
    transformed_gtf: ./ref_files/hg38/hg38_gtf_full_table.txt.gz
    GPL: ./ref_files/GPL/
  ```

  download from:  
  - gencode V29: <https://www.gencodegenes.org/human/release_29.html>  
  - ENOCDE: <https://www.encodeproject.org/references/ENCSR151GDH/>
  - STAR index: <https://www.encodeproject.org/files/ENCFF598IDH/>
  - rsem index: <https://www.encodeproject.org/files/ENCFF285DRD/>
  - gtf: <https://www.encodeproject.org/files/gencode.v29.primary_assembly.annotation_UCSC_names/>  
  - genome_fasta:  <https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/>
  - transcriptome_fasta: <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz>
  - chromosome_size: <https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/>

- mm10:  
  
  ```yaml
  mm10:
    salmon_index: ./ref_files/mm10/salmon_index
    STAR_index: ./ref_files/mm10/STAR_index
    gtf: ./ref_files/mm10/gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz
    transformed_gtf: ./ref_files/mm10/mm10_gtf_full_table.txt.gz
    GPL: ./ref_files/GPL/
  ```

  download from:  
  - gencode M21: <https://www.gencodegenes.org/mouse/release_M21.html>  
  - ENCODE: <https://www.encodeproject.org/references/ENCSR496QMW/>
  - STAR index: <https://www.encodeproject.org/files/ENCFF795ZFF/>
  - rsem index: <https://www.encodeproject.org/files/ENCFF363TFV/>
  - gtf: <https://www.encodeproject.org/files/gencode.vM21.primary_assembly.annotation_UCSC_names/>  
  - genome_fasta: <https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/>  
  - transcriptome_fasta: <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.transcripts.fa.gz>  
  - chromosome_size: <https://www.encodeproject.org/files/mm10_no_alt.chrom.sizes/>  
