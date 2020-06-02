#!/bin/bash

# UMI=10
# CBC=12
# CELLS=5000
# TOPCELLS=10000
# METHOD=directional
# HAMMING=1
# UMI_THRESH=30
# OVERLAP_THRESH=10
# READLEN=25
# MISMATCHES=1
# THRESH=0.005

DIR="/seq/lincRNA/Projects/MIP/mip-pipeline"

source /broad/software/scripts/useuse
reuse .python-3.5.1
source /seq/lincRNA/Ben/VENV_MIP/bin/activate
reuse Bowtie
reuse Samtools
reuse FASTX-Toolkit
reuse BEDTools

SAMPLESHEET=$1
CONFIG=$2 
# example config: /seq/lincRNA/Projects/MIP/mip-pipeline/file_input/config.json

mkdir -p logs

# snakemake -s "$DIR/file_input/Snakefile" -j 100 -r -p --cluster "$DIR/bin/snake-qsub" -k --config sample_sheet=$SAMPLESHEET umi_thresh=$UMI_THRESH overlap_thresh=$OVERLAP_THRESH umi=$UMI cbc=$CBC cells=$CELLS topcells=$TOPCELLS method=$METHOD hamming=$HAMMING readlen=$READLEN mismatches=$MISMATCHES -w 60
# snakemake -s "$DIR/bin/Snakefile" -j 100 -r -p --cluster "source /broad/software/scripts/useuse ; reuse -q GridEngine8 ; qsub -q gsa -j y -cwd -o qsub_output -e qsub_error -l virtual_free={cluster.mem}G" --jobscript "/seq/lincRNA/Projects/MIP/mip-pipeline/file_input/jobscript.sh" --cluster-config "/seq/lincRNA/Projects/MIP/mip-pipeline/file_input/cluster.json" -k --config sample_sheet=$SAMPLESHEET umi_thresh=$UMI_THRESH overlap_thresh=$OVERLAP_THRESH umi=$UMI cbc=$CBC method=$METHOD hamming=$HAMMING readlen=$READLEN mismatches=$MISMATCHES thresh=$THRESH -w 60 # cells=$CELLS topcells=$TOPCELLS
snakemake -s "$DIR/bin/Snakefile" -j 100 -r -p --cluster "source /broad/software/scripts/useuse ; reuse -q GridEngine8 ; qsub -q gsa -j y -cwd -o {cluster.output} -l virtual_free={cluster.memory}G" --jobscript "/seq/lincRNA/Projects/MIP/mip-pipeline/file_input/jobscript.sh" --cluster-config "/seq/lincRNA/Projects/MIP/mip-pipeline/file_input/cluster.json" -k --configfile $CONFIG --config sample_sheet=$SAMPLESHEET -w 60
