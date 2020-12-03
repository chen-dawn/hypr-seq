# HyPR-seq

This repo contains code to analyze HyPR-seq data, an assay introduced by Marshall et al. 2020 (https://www.biorxiv.org/content/10.1101/2020.06.01.128314v1). It takes fastqs and returns count tables, as well as some stats and QC plots along the way. Code for designing HyPR-seq probes is also included, as well as code to process some of the data we analyzed in the paper, although this will be fleshed out more later.

## Analyzing HyPR-seq data
To summarize, we first filter reads, build a whitelist of real cells, extract cell barcodes and UMIs, map our reads, collapsed bead barcode doublets, and return a count matrix. We use UMI-Tools for dealing with UMIs, fastp for filtering, bowtie for mapping, and Snakemake for pipelining. Currently, the code runs off a Snakefile which is configured to submit jobs to an SGE cluster (although we are trying to move this to gcp!). 

Code for most of the analyses is in bin/, with cluster specific files in cluster-config.

#### Snakefile for analyzing HyPR-seq data
- Install miniconda and the following conda environment:  conda env create --file cluster-config/hyprseq.yml 
- Copy or edit the cluster-config/jobscript.sh file as needed for your cluster environment
- Set up the Samplesheet, Probe list, and config files as described below (see examples in example/)
- Edit/run the following commands (example uses SLURM cluster):

  git clone --branch Sherlock git@github.com:EngreitzLab/hypr-seq.git
  
  snakemake -s Snakefile -j <MaxJobs> -r -p \
    --directory <YourProjectDirectory> \
    --cluster "sbatch -J HYPR -n 1 -c 1 --export=ALL --mem {cluster.memory}G --wrap" \
    --jobscript cluster-config/jobscript.sh \
    --cluster-config cluster-config/cluster.json \
    -k -w 60 --configfile example/params.json --config sample_sheet=<SampleSheet> codedir=bin/

#### Samplesheet specs
See example in example/SampleList.txt

#### Config specs
See example in example/params.json

#### Probe List
See example in example/ProbeListExample.txt

## Designing HyPR-seq probes
Code for designing HyPR-seq probes is in probe-design. 

## Paper-specific analyses
Code in paper_analyses (currently just files to process the CRISPR screen and kidney experiments) is used to make figures that we used in our paper. These will be updated to be more fleshed out.
