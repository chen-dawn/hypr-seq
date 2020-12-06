# Paper Analyses

This repo contains code to reproduce the trickier analyses from the paper (in particular, the kidney and splenocyte single-cell datasets) and to analyze HyPR-seq CRISPR screens.

## Kidney
kidney_pipeline.py has code to cluster the HyPR-seq data, assign cell types, and plot UMAPs, heatmaps, and cell-type frequency changes.

## Spleen
splenocyte_pipeline.py has code to compare the 250- and 1000-probe experiments, assign HyPR-seq cells to clusters based on pre-existing scRNA-seq data, and to plot gene expression heatmaps.

## CRISPR Screen
analyze_crispr_screen.py has code used to process HyPR-seq data from the small pooled CRISPR screen performed in the paper. It can assign guides to cells, aggregate cells by guide, and plot relative gene expression. To run, it requires the following files:
- Collapsed count table, {sample}/{sample}\_collapsed\_emul.count\_matrix.txt from the output folder
- Uncollapsed count table, {sample}/{sample}\_uncollapsed\_emul.count\_matrix.txt from the output folder
- Barcode stats file, {sample}/tmp/{sample}.collapse\_bc\_stats.txt from the output folder
- Guide info file, tab-separated file with guide info (guide name, associated barcode, whether it is a negative control, etc.) (see example)
- Target genes, comma-separated list of genes to plot (i.e. 'GATA1' or 'GATA1,HDAC6,BFP')
- Housekeeping genes, comma-separated list of genes to use as normalizing factors (i.e. 'GAPDH' or 'GAPDH,RPL13A,ACTB')
- Outdir, path to directory to output the results
- Name, prefix for the output files
- Codedir, path to cloned hypr-seq directory

Optional
- Chastity, how strict to be when identifying expressed barcodes in cells (default: 6)
- Entropy threshold, the maximum "entropy" (of the distribution of barcode probe counts per cell) that counts as a single cell, purely historical (default: 0.5)
- dCas9 filter, if included, will filter out cells with low (bottom 15%) expression of KRAB-dCas9 (default: False)
- dCas9 transcript, what probe name to use as the proxy for KRAB-dCas9 expression (default: 'BFP')


python analyze_crispr_screen.py --collapsed \<collapsed_file\> --uncollapsed \<uncollapsed_file\> --barcode_stats \<barcode_stats_file\> \
	--guide_info \<collapsed_file\> --target_genes \<GENE1,GENE2\> --housekeeping_genes \<GENE1,GENE2,GENE3\> --outdir \<OutputDirectory\> \
	--name \<Sample\> --codedir \</path/to/hypr-seq\>
