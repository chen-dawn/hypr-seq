# HyPR-seq

This repo contains code to analyze HyPR-seq data, an assay introduced by Marshall et al. 2020 (https://www.biorxiv.org/content/10.1101/2020.06.01.128314v1). It takes fastqs and returns count tables, as well as some stats and QC plots along the way. Code for designing HyPR-seq probes is also included, as well as code to process some of the data we analyzed in the paper, although this will be fleshed out more later.

## Analyzing HyPR-seq data
To summarize, we first filter reads, build a whitelist of real cells, extract cell barcodes and UMIs, map our reads, collapsed bead barcode doublets, and return a count matrix. We use UMI-Tools for dealing with UMIs, fastp for filtering, bowtie for mapping, and Snakemake for pipelining. Currently, the code runs off a Snakefile which is configured to submit jobs to an SGE cluster (although we are trying to move this to gcp!). 

Code for most of the analyses is in bin/, with cluster specific files in cluster-config.

#### Samplesheet specs

#### Config specs

## Designing HyPR-seq probes
Code for designing HyPR-seq probes is in probe-design. 

## Paper-specific analyses
Code in paper_analyses (currently just files to process the CRISPR screen and kidney experiments) is used to make figures that we used in our paper. These will be updated to be more fleshed out.
