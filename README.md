# HyPR-seq

This repo contains code to analyze HyPR-seq data, an assay introduced by Marshall et al. 2020 (https://www.biorxiv.org/content/10.1101/2020.06.01.128314v1). It takes fastqs and returns count tables, as well as some stats and QC plots along the way. Code for designing HyPR-seq probes is also included, as well as code to process some of the data we analyzed in the paper, although this will be fleshed out more later.

## Analyzing HyPR-seq data
To summarize, we first filter reads, build a whitelist of real cells, extract cell barcodes and UMIs, map our reads, collapsed bead barcode doublets, and return a count matrix. We use UMI-Tools for dealing with UMIs, fastp for filtering, bowtie for mapping, and Snakemake for pipelining. We have configured the code to submit jobs to a cluster scheduler (currently we have gotten it to work on SGE and SLURM, and we are trying to move to gcp).  

Code for going from fastqs to count tables (+ QCs) is in bin/, cluster + conda files in cluster-config/. The files in paper_analyses/ contain code to make the figures in the paper (the kidney + splenocyte experiments), as well as scripts for analyzing the CRISPR screen (if you are interested in conducting a screen of your own, check that directory and the README therein). The example/ directory contains sample files for the input (see below for format). If you are looking to design probes to genes of interest, see probe-design/ and the README.

#### Snakefile for analyzing HyPR-seq data
- Install miniconda and the following conda environment:  conda env create --file cluster-config/hyprseq.yml 
- Copy or edit the cluster-config/jobscript.sh file as needed for your cluster environment
- Set up the Samplesheet, Probe list, and config files as described below (see examples in example/)
- Edit/run the following commands (example uses SLURM cluster):

  git clone --branch Sherlock git@github.com:EngreitzLab/hypr-seq.git
  
  snakemake -s Snakefile -j \<MaxJobs\> -r -p \
    --directory <YourProjectDirectory> \
    --cluster "sbatch -J HYPR -n 1 -c 1 --export=ALL --mem {cluster.memory}G --wrap" \
    --jobscript cluster-config/jobscript.sh \
    --cluster-config cluster-config/cluster.json \
    -k -w 60 --configfile example/params.json --config sample_sheet=<SampleSheet> codedir=bin/

(Note, if you are interested in running on an SGE cluster, replace the --cluster flag with something like "qsub -q gsa -j y -cwd -o {cluster.output} -l virtual_free={cluster.memory}G".)

#### Samplesheet specs
See example in example/SampleList.txt

The samplesheet is a tab-delimited text file, with one line per sample (pair of FASTQ files), with the following columns:
- ID - **unique** integer ID per sample, must not be duplicated on one sample sheet (so you can have the same English name for different samples, if you want)
- Name - plain English name describing the sample (no spaces), can be shared between samples (outputs will be stored in a directory called Name-ID)
- Read1 - full path to the fastq.gz file containing read1 (UMI + transcript-binding probe)
- Read2 - full path to the fastq.gz file containing read2 (cell barcode)
- ProbeName - plain name to describe the bowtie index that will be built from the provided probe list (helpful to name it in agreement with ProbePath below, but it doesn't matter)
- ProbePath - full path to the probe list file (see specs below)
- CellNum - how many cells to consider in the analysis, either an integer or the string "Guess" (in which case UMI-tools will use the knee plot to infer the number of cells); we recommend starting with "Guess" and if it's wrong based on manual inspection re-running with your best guess
- ReadsToUse [OPTIONAL] - how many reads to include in the analysis (to facilitate downsampling) if you want all the reads in the file, enter 0 (or omit the column)

#### Config specs
See example in example/params.json

The config file is a JSON file containing parameters used for the analysis. Most of these parameters are baked into the final version of the library design and will not need to be changed, except potentially "thresh" (yes, bad name, I know). The required parameters (and their default values in example/) are:
- umi: UMI length (default: 10)
- cbc: cell-barcode length (default: 12)
- readlen: probe binding site length (default: 25)
- mismatches: the number of mismatches for read mapping using Bowtie (default: 1)
- hamming: Hamming distance threshold for UMI deduping in UMI-tools (default: 1)
- method: UMI-tools specific dedup method (default: 'directional')
- thresh: fractional threshold for identifying pairs of cell barcodes that may have originated in the same droplet (see below) (default: 0.005)

The "thresh" parameter is the only parameter we adjust on a semi-regular basis. For more information on our cell-barcode deduping procedure, see below. We recommend starting with the default value and increasing it as needed.

#### Probe List
See example in example/ProbeListExample.txt

The probe list file is a tab-delimited file with two columns, the probe name and probe sequence. **The probe names should be in the form [gene_name]5P_[probe_number]**, and the probe sequence should be the 25bp unique RNA-binding sequence of the 5' probe. 

## Designing HyPR-seq probes
Code for designing HyPR-seq probes is in probe-design. 

## Paper-specific analyses
Code in paper_analyses (currently just files to process the CRISPR screen and kidney/splenocyte experiments) is used to make figures that we used in our paper. These will be updated shortly to be more fleshed out!
