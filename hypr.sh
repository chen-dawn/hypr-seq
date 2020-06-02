#!/bin/bash

DIR="/seq/lincRNA/Ben/bin/hypr-seq"

source /broad/software/scripts/useuse
reuse .python-3.5.1
source /seq/lincRNA/Ben/VENV_MIP/bin/activate
reuse Bowtie
reuse Samtools
reuse FASTX-Toolkit
reuse BEDTools

SAMPLESHEET=$1
CONFIG=$2 

mkdir -p logs

snakemake -s "$DIR/Snakefile" -j 100 -r -p --cluster "source /broad/software/scripts/useuse ; reuse -q GridEngine8 ; qsub -q gsa -j y -cwd -o {cluster.output} -l virtual_free={cluster.memory}G" --jobscript "/seq/lincRNA/Projects/MIP/mip-pipeline/file_input/jobscript.sh" --cluster-config "/seq/lincRNA/Projects/MIP/mip-pipeline/file_input/cluster.json" -k --configfile $CONFIG --config sample_sheet=$SAMPLESHEET -w 60
