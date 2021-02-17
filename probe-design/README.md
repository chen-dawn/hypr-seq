# Probe Design

This repo contains code to design HyPR-seq probes

## Requirements
In addition to everything in hyprseq.yml (in cluster-config), the probe designer (for historical reasons) also relies on Matlab and an older version of BLAST. We are interested in re-writing the code to not rely on these, but for now, the code needs Matlab 2017a (I think newer ones will also work, potentially up to 2020a, but I've only tested on 2017a and the function "blastlocal" still has to exist) and BLAST 2.2.17 (which is currently not supported but can be found here: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.17/). We are working on changing this.

## Config file
An example config file can be found in this directory. Basically you need to specify:
- input_table: two column file with the gene names you want to design probes for repeated in both columns (historical, we can probably change this)
- gtf: gtf file for identifying coordinates of genes of interest (NOTE: different gtf files I have used have different formatting (in particular, my ones for mouse and human had different tags and such), I believe this code should be extensible to the formats I've seen, but if you have differently formatted gtfs you might have to edit. ping me at bgrd [at] stanford [dot] edu for an example gtf or submit a pull request and we can improve the gtf parsing)
- genome_fasta: path to genome fasta file for extracting sequence (also need .fai in same directory)
- intons: boolean, do you want to design intron probes as well
- blastdb: past to a BLAST db for running local blast searches (should end in .fna, although there are associated other files that are necessary too, I think), can be downloaded from ncbi (e.g. ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_rna.fna.gz)
- probes_per_gene: path to a probe_number file, simply a two column tab delimited file with no header where the first column is genes and the second is how many probes to design (optional, if excluded or gene not in, defaults to 6 probes)
- gene_name_field: where is the gene name stored in the gtf file (gene_id, gene_name, etc.) (optional, defaults to "gene_id")
- search_tag_basic: conditional gtf file parsing option, sometimes a "basic" tag marks the canonical isoform, do you want to limit your search to this (boolean) (optional, defaults to false)

## Example
Included here is a log file for generating probes by calling the snakemake. 

## Outputs
Primarily, a table called probes_for_idt.txt, which contains a table that can be directly uploaded to IDT for ordering. If you make any modifications to the molecular biology and need to adjust the flanking sequences, you will have to dig into HCRProbeGeneration.m (it should be straightforward what is being appended). A list of how many successful probes per gene is in probe_stats.txt, and a list of genes that failed is in genes_with_no_probes.txt (probes.csv is an intermediate and can be deleted). It's always a good idea to blat some of the probe binding sites as a sanity check before placing a big order.
