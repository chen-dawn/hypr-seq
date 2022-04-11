#!/bin/sh
# properties = {properties}

## Add any needed cluster setup commands here
source /broad/software/scripts/useuse
# source $HOME/.bioinfo
use .bowtie2-2.3.4.3
use .bamtools-2.5.1
use .bedtools-2.29.0

use Anaconda3
source activate hyprseq

{exec_job}
