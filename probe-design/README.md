# Probe Design

This repo contains code to design HyPR-seq probes

## Requirements
In addition to everything in hyprseq.yml (in cluster-config), the probe designer (for historical reasons) also relies on Matlab and an older version of BLAST. We are interested in re-writing the code to not rely on these, but for now, the code needs Matlab 2017a (newer ones might work, but not all new ones, the function "blastlocal" still has to exist) and BLAST 2.2.17 (which is currently not supported but can be found here: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.17/). We are working on changing this.

## Config file

## Example

