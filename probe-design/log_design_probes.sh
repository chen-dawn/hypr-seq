#!/bin/bash

# NEED TO LOAD MATLAB AND BLAST
ml matlab/R2017a
export PATH=$PATH:/home/users/bgrd/bin/blast-2.2.17/bin/

# activate conda environment
conda activate hyprseq

CONFIG=210203_probes.json

mkdir -p logs HomologyRegionsHCR

# run snakemake
snakemake -s hypr-seq/probe-design/Snakefile -j 10 -r -p \
--directory <MY_DIR> \
--cluster "sbatch -J HYPR -n 1 -c 1 --export=ALL --mem {cluster.memory}G --wrap" \
--jobscript hypr-seq/cluster-config/jobscript.sh \
--cluster-config hypr-seq/cluster-config/cluster.json \
-k -w 60 --configfile $CONFIG --config codedir=hypr-seq/probe-design/ 