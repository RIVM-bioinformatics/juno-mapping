[![Super linter](https://github.com/RIVM-bioinformatics/juno-mapping/actions/workflows/super-linter.yaml/badge.svg)](https://github.com/RIVM-bioinformatics/juno-mapping/actions/workflows/super-linter.yaml)


<div align="center">
    <h1>Juno-Mapping</h1>
    <br />
    <h2>Pipeline for reference-based variant calling of bacterial genomes.</h2>
    <br />
</div>

## Pipeline information

* **Authors:**            Boas van der Putten.
* **Organization:**         National Institute for Public Health and the Environment (RIVM)
* **Department:**           Centre for Research Infectious Diseases Diagnostics and Screening (IDS), Informatie en Beheer (IBR)
* **Start date:**           01 - 12 - 2022
* **Commissioned by:**      Richard Anthony

## About this project

The goal of this pipeline is to generate variant calls from raw fastq files. The input of the pipeline is raw Illumina paired-end data (read length > 99 bp) in the form of two fastq files (with extension .fastq, .fastq.gz, .fq or .fq.gz), containing the forward and the reversed reads ('R1' and 'R2' must be part of the file name, respectively). On the basis of the generated quality reports, low quality and contaminated samples can be excluded for downstream analysis. __Note:__ The pipeline is developed mainly with _Mycobacterium tuberculosis_ in mind, but could theoretically be applied to other species.

The pipeline uses the following tools:
1. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (Andrews, 2010) is used to assess the quality of the raw Illumina reads
2. [FastP](https://github.com/OpenGene/fastp) (Chen, Zhou, Chen and Gu, 2018) is used to remove poor quality data and adapter sequences 
3. [Picard](https://broadinstitute.github.io/picard/) determines the library fragment lengths
4. [bwa](https://bio-bwa.sourceforge.net/bwa.shtml) for read mapping
5. [GATK](https://gatk.broadinstitute.org/hc/en-us) for variant calling (Mutect2) and various VCF analyses.
6. [bcftools](https://samtools.github.io/bcftools/bcftools.html) for VCF filtering.
7. [samtools](https://github.com/samtools/samtools) for processing of SAM and BAM files.
8. [Bbtools](https://jgi.doe.gov/data-and-tools/bbtools/) (Bushnell, 2014) is used to generate scaffold alignment metrics 
9. [MultiQC](https://multiqc.info/) (Ewels et al., 2016) is used to summarize analysis results and quality assessments in a single report for dynamic visualization.
10. [Kraken2](https://ccb.jhu.edu/software/kraken2/) (Wood, Lu and Langmead, 2019) and [Bracken](http://ccb.jhu.edu/software/bracken/) (Lu, Breitwieser, Thielen and Salzberg, 2016) for identification of bacterial species.  

![Image of pipeline](https://github.com/RIVM-bioinformatics/juno-assembly/blob/master/files/juno_mapping.svg)

## Prerequisities

* **Linux + conda/mamba** A Linux-like environment with at least 'miniconda' installed. See below for instructions on how to install mamba through conda.
* Preferentially **[Singularity](https://sylabs.io/guides/latest/user-guide/)** installed and working. If you are not using singularity, you have to run the pipeline using the argument `--no-containers`.


## Installation

**IMPORTANT NOTE**: You need to have [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/en/latest/) installed and, preferably, also [Singularity](https://sylabs.io/guides/latest/user-guide/) so that every step is containerized and therefore more reproducible. There is an option to run without using singularity (option `--no-containers`) but conda or mamba are mandatory. 

1. Clone the repository:

```bash
git clone https://github.com/RIVM-bioinformatics/juno-mapping.git
```
Alternatively, you can download it manually as a zip file (you will need to unzip it then). If you decide to download the zip only, the pipeline version will not be stored in the audit trail.

2. Enter the directory with the pipeline and install the pipeline. The script will check whether mamba is available and installs it from mamba.yaml if unavailable:

```bash
cd juno-mapping
bash install_juno_mapping.sh
```

## Parameters & Usage

### Command for help

* ```-h, --help``` Shows the help of the pipeline:

```
usage: juno_mapping.py [-h] -i DIR [-o DIR] [-w DIR] [-ex FILE] [-p PATH] [-l] [-tl INT] [-u] [-n] [-q QUEUE] [--no-containers] [--snakemake-args [SNAKEMAKE_ARGS ...]] -s STR [--reference FILE] [--mask FILE]
                       [--disable-mask] [--db-dir DIR] [-mpt INT] [-ws INT] [-ml INT] [-md INT] [-maf FLOAT]

Juno-mapping pipeline for reference mapping of bacterial genomes

options:
  -h, --help            show this help message and exit
  -i DIR, --input DIR   Relative or absolute path to the input directory. It must contain all the raw reads (fastq) files for all samples to be processed (not in subfolders).
  -o DIR, --output DIR  Relative or absolute path to the output directory. If none is given, an 'output' directory will be created in the current directory.
  -w DIR, --workdir DIR
                        Relative or absolute path to the working directory. If none is given, the current directory is used.
  -ex FILE, --exclusionfile FILE
                        Path to the file that contains samplenames to be excluded.
  -p PATH, --prefix PATH
                        Conda or singularity prefix. Basically a path to the place where you want to store the conda environments or the singularity images.
  -l, --local           If this flag is present, the pipeline will be run locally (not attempting to send the jobs to an HPC cluster**). The default is to assume that you are working on a cluster. **Note that currently
                        only LSF clusters are supported.
  -tl INT, --time-limit INT
                        Time limit per job in minutes (passed as -W argument to bsub). Jobs will be killed if not finished in this time.
  -u, --unlock          Unlock output directory (passed to snakemake).
  -n, --dryrun          Dry run printing steps to be taken in the pipeline without actually running it (passed to snakemake).
  -q QUEUE, --queue QUEUE
                        Name of the queue that the job will be submitted to if working on a cluster.
  --no-containers       Use conda environments instead of containers.
  --snakemake-args [SNAKEMAKE_ARGS ...]
                        Extra arguments to be passed to snakemake API (https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html).
  -s STR, --species STR
                        Species to use, choose from: ['mycobacterium_tuberculosis']
  --reference FILE      Reference genome to use. Default is chosen based on species argument, defaults per species can be found in: /mnt/db/juno/mapping/[species]
  --mask FILE           Mask file to use, defaults per species can be found in: /mnt/db/juno/mapping/[species]
  --disable-mask        Disable masking, use at your own risk: this might cause the appearance of low quality variants in the final VCF file.
  --db-dir DIR          Kraken2 database directory.
  -mpt INT, --mean-quality-threshold INT
                        Phred score to be used as threshold for cleaning (filtering) fastq files.
  -ws INT, --window-size INT
                        Window size to use for cleaning (filtering) fastq files.
  -ml INT, --minimum-length INT
                        Minimum length for fastq reads to be kept after trimming.
  -md INT, --minimum-depth-variant INT
                        Minimum length for fastq reads to be kept after trimming.
  -maf FLOAT, --minimum-allele-frequency FLOAT
                        Minimum allele frequency to filter variants on.
```

### The base command to run this program. 

Before running the pipeline you should always activate the conda environment.

```bash
conda activate juno_mapping
python juno_mapping.py -i [path/to/input/dir] -s [species] -o [path/to/output/dir]
```

You can deactivate the environment once you are finished.

```bash
conda deactivate
```

### An example on how to run the pipeline.

```bash
conda activate juno_mapping
python juno_mapping.py -i my_data -o my_results -s Mycobacterium_tuberculosis --local
conda deactivate
```

## Explanation of the output

* `log`: Log files with output and error files from each Snakemake rule/step that is performed. 
* `audit_trail`: Information about the versions of software and databases used.
* **output per sample:** The pipeline will create one subfolder per each step performed. These subfolders will in turn contain the results per sample and/or one with results combined for all samples (if applicable). To understand the ouptut, please refer to the manuals of each individual tool. The generated subfolders are:
    - clean_fastq: 
    - de_novo_assembly: results from the assembly step without filtering out the small contigs.
    - de_novo_assembly_filtered: results from the assembly step containing only contigs larger than 500bp. This is what you would normally use in downstream analyses that require assembly (fasta files) as input.
    - multiqc: 
    - qc_clean_fastq: 
    - qc_de_novo_assembly: results of the quality control of the assemblies. Includes results of different tools (refer to those tools for interpretation).
    - qc_raw_fastq:  
    - identify_species:   

    - `clean_fastq`: contains the fastq files after filtering low quality reads and trimming ends and/or adapters. This is what you would use in downstream analyses that require fastq files as input.
    - `identify_species`: results of bracken for species identification. The main result is a file called top1_species_multireport.csv which summarizes the main species found in the sample.
    - `mapped_reads`: Mapped reads in BAM format, with indices.
    - `multiqc`: MultiQC report for all samples run.
    - `qc_clean_fastq`: results of the quality control for fastq files. Run _after_ filtering and trimming.
    - `qc_mapping`: results of the quality control for mapped BAM files.
    - `qc_raw_fastq`: results of the quality control for fastq files. Run _before_ filtering and trimming.
    - `qc_variant_calling`: results of the quality control for VCF files.
    - `reference`: the reference used with corresponding indices.
    - `variants`: filtered VCF files, comprisong the final output of this pipeline.
    - `variants_raw`: several versions of unfiltered VCF files. Useful for quality control or if you're interested in low abundance variants.
    - `variants_snps_only`: filtered VCF files containing only SNPs.

## Issues  

* The pipeline currently only supports LSF clusters or local running.
* Any issue can be reported in the [Issues section](https://github.com/RIVM-bioinformatics/juno-mapping/issues) of this repository.

## License
This pipeline is licensed with an AGPL3 license. Detailed information can be found inside the 'LICENSE' file in this repository.

## Contact
* **Contact person:**       Boas van der Putten
* **Email**                 ids-bioinformatics@rivm.nl


## Contribution guidelines
Juno pipelines use a [feature branch workflow](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow). To work on features, create a branch from the `main` branch to make changes to. This branch can be merged to the main branch via a pull request. Hotfixes for bugs can be committed to the `main` branch.

Please adhere to the [conventional commits](https://www.conventionalcommits.org/) specification for commit messages. These commit messages can be picked up by [release please](https://github.com/googleapis/release-please) to create meaningful release messages.

This pipeline is based on Juno-template.