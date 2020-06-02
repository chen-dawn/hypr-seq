#!/bin/bash

## USAGE
# run_probe_design.sh <config.json> 

# Config:
# GeneList.xlsx
# Genome fasta
# Gene .gtf
# Introns T/F

DIR=/seq/lincRNA/Projects/MIP/mip-pipeline/probe_design/

source /broad/software/scripts/useuse
reuse .python-3.5.1
source /seq/lincRNA/Ben/VENV_MIP/bin/activate

CONFIG=$1

mkdir -p logs HomologyRegionsHCR

snakemake -s "$DIR/Snakefile" -j 100 -r -p --cluster "source /broad/software/scripts/useuse ; reuse -q GridEngine8 ; qsub -q gsa -j y -cwd -o {cluster.output} -l virtual_free={cluster.memory}G" --jobscript "/seq/lincRNA/Projects/MIP/mip-pipeline/file_input/jobscript.sh" --cluster-config "/seq/lincRNA/Projects/MIP/mip-pipeline/file_input/cluster.json" -k --configfile $CONFIG -w 60

