#!/usr/bin/env python

## USAGE python analyze_count.py <count_file> <count_plots> <count_log> <whitelist>
# takes in the output from umi-tools count command and does things with plotting various mips

import sys
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')
import itertools
import distance
import argparse
from scipy.stats import zscore
import openpyxl

import mip_tools

maxumispercell = 20000
minumispercell = 0

def analyze_count(flat_counts_file, plots, log, stats, 
                  collapsed, barcodes_stats_file,
                  cell_x_gene_file, cell_x_gene_rescaled_file, cell_x_probe_file, cell_x_probe_rescaled_file,
                  sum_probes_file, sum_genes_file, average_probes_file, average_genes_file,
                  zscores_gene_file, zscores_probe_file, corr_file, problematic_probes_file,
                  excel_counts_file, excel_probes_file, excel_gene_counts_file, excel_zscores_file):
    
    # read flat file
    raw_flat = pd.read_table(flat_counts_file)

    # pivot to genes x cells and then transpose
    raw_probes = raw_flat.pivot(index='gene', columns='cell', values='count').fillna(0).T

    # filter out weird "cells"
    filtered_probes = mip_tools.filter_cells(raw_probes)

    # check whether this is collapsed or uncollapsed file and if it is remove giant droplets
    if barcodes_stats_file:
        bc_stats = pd.read_table(barcodes_stats_file, index_col=0)
        bcs_to_keep = bc_stats.loc[(bc_stats.Singleton == True) | (bc_stats.Multiplet == True)].index
        filtered_probes = filtered_probes.loc[filtered_probes.index.intersection(bcs_to_keep)]

    # stats_suffix = '_collapsed' if collapsed else '_uncollapsed'
    stats_suffix = '_'+collapsed


    stats.write('Cells_after_filtering{}\t{}\n'.format(stats_suffix, filtered_probes.shape[0]))
    stats.write('Probes{}\t{}\n'.format(stats_suffix, filtered_probes.shape[1]))
    stats.write('Dropped_Cells{}\t{}\n'.format(stats_suffix, len(raw_probes) - len(filtered_probes)))

    # get in order
    filtered_probes = filtered_probes.reindex(sorted(filtered_probes.columns), axis=1)

    # get to constant UMIs/cell
    filtered_probes_rescaled = mip_tools.rescale_matrix_to_constant_row_sums(filtered_probes)

    # write counts per probe and rescaled
    filtered_probes.to_csv(cell_x_probe_file, sep='\t', header=True, index=True)
    filtered_probes_rescaled.to_csv(cell_x_probe_rescaled_file, sep='\t', header=True, index=True)

    # compute and plot sum
    mip_tools.sum_across_cells_and_plot(filtered_probes, plots, sum_probes_file, None, title='Probes')

    # compute and plot average
    mip_tools.average_across_cells_and_plot(filtered_probes_rescaled, plots, average_probes_file, title='Probes')

    # compute correlation (on rescaled matrix)
    probe_corr = mip_tools.compute_and_plot_correlation(filtered_probes_rescaled, plots, corr_file)

    # average correlation within genes and compare to abundances
    mip_tools.average_correlation_and_compare_to_expression(probe_corr, filtered_probes_rescaled, plots)

    # compute z scores on probes
    mip_tools.calculate_zscores(filtered_probes_rescaled, zscores_probe_file, None, log)

    # write excel file with the counts per probe
    mip_tools.write_counts_per_probe_excel(filtered_probes, excel_probes_file)

    # highlight problematic probes
    # would be very easy to build in the correlation stuff here, although I won't
    mip_tools.find_problematic_probes(filtered_probes, plots, stats, stats_suffix, problematic_probes_file)

    # collapse individual probes to transcripts
    filtered_genes = pd.DataFrame(index=filtered_probes.index)

    all_probes = filtered_probes.columns

    probe_targets = sorted(list({mip_tools.extract_base_probe_name(probe) for probe in all_probes}))

    # collapse probes by finding all columns that match regex and sum across
    # might also want to swap this out for the sum with exclusion list function in mip_tools
    for probe in probe_targets:
        if probe in mip_tools.mip_barcodes:
            # needs regex
            filtered_genes[probe] = filtered_probes.filter(regex='^{}5P_[0-9]+'.format(probe), axis=1).sum(axis=1)
        else:
            # treat normally
            filtered_genes[probe] = mip_tools.get_probes_for_gene(filtered_probes, probe).sum(axis=1)
        
    stats.write('Genes{}\t{}\n'.format(stats_suffix, filtered_genes.shape[1]))
    
    # reindex column names
    filtered_genes = filtered_genes.reindex(sorted(filtered_genes.columns), axis=1)

    # get to constant UMIs/cell
    filtered_genes_rescaled = mip_tools.rescale_matrix_to_constant_row_sums(filtered_genes)

    # write counts per gene and rescaled
    filtered_genes.to_csv(cell_x_gene_file, sep='\t', header=True, index=True)
    filtered_genes_rescaled.to_csv(cell_x_gene_rescaled_file, sep='\t', header=True, index=True)

    # compute and plot sum
    mip_tools.sum_across_cells_and_plot(filtered_genes, plots, sum_genes_file, excel_gene_counts_file, title='Genes')

    # compute and plot average
    mip_tools.average_across_cells_and_plot(filtered_genes_rescaled, plots, average_genes_file, title='Genes')

    # compute z scores on genes
    mip_tools.calculate_zscores(filtered_genes_rescaled, zscores_gene_file, excel_zscores_file, log)

    # # plot and save counts per gene
    # mip_tools.plot_counts_per_gene(collapsed_counts, plots, gene_counts, excel_gene_counts_file)

    # plots per gene
    # do this for both raw counts and rescaled counts
    for probe in probe_targets:
        # plot histogram and cdf for each probe
        mip_tools.plot_probe(filtered_genes, probe, plots)
        mip_tools.plot_probe(filtered_genes_rescaled, probe, plots, title='Rescaled')

        # also plot the different probes comprising each probe (don't need to do this for the rescaled ones)
        mip_tools.plot_probes_per_transcript(filtered_probes, probe, plots)

        # also do pairplot for all the probes (to show correlation)
        mip_tools.plot_pairplot_of_individual_probes(filtered_probes, probe, plots)
        mip_tools.plot_pairplot_of_individual_probes(filtered_probes_rescaled, probe, plots, title='Rescaled')

    # plot umis/cell
    umis_per_cell = filtered_genes.sum(axis=1).values
    stats.write('Avg_UMIs/cell_after_filtering{}\t{}\n'.format(stats_suffix, np.mean(umis_per_cell)))
    stats.write('Med_UMIs/cell_after_filtering{}\t{}\n'.format(stats_suffix, np.median(umis_per_cell)))

    mip_tools.plot_umis_per_cell_analyze_count(umis_per_cell, plots)

    # write to excel files for Jamie
    try:
        filtered_genes.T.reset_index().to_excel(excel_counts_file, header=False, index=False)
    except:
        log.write("Too many CBCs, can't write excel file in this direction, try transposed\n")
        filtered_genes.to_excel(excel_counts_file)

    log.write('Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze the count matrix from a Hypr-seq experiment and return useful plots.')
    
    parser.add_argument("--flat", dest="flat_file", type=str,  help="Raw counts file (flat and not filtered)")
    parser.add_argument("--plots", dest="plots_file", type=str,  help="Where to output plots")
    parser.add_argument("--log", dest="log_file", type=str,  help="Where to output log")
    parser.add_argument("--stats", dest="stats_file", type=str,  help="Where to output stats")
    parser.add_argument("--cell_x_gene", dest="cell_x_gene_file", type=str,  help="Count table (cells x genes)")
    parser.add_argument("--cell_x_gene_rescaled", dest="cell_x_gene_rescaled_file", type=str,  help="Count table (cells x genes) normalized row sums")
    parser.add_argument("--cell_x_probe", dest="cell_x_probe_file", type=str,  help="Count table (cells x probes)")
    parser.add_argument("--cell_x_probe_rescaled", dest="cell_x_probe_rescaled_file", type=str,  help="Count table (cells x probes) normalized row sums")
    parser.add_argument("--sum_probes", dest="sum_probes_file", type=str, help="Sum of probes across cells")
    parser.add_argument("--sum_genes", dest="sum_genes_file", type=str, help="Sum of genes across cells")
    parser.add_argument("--average_probes", dest="average_probes_file", type=str, help="Average of probes across cells (after scaling)")
    parser.add_argument("--average_genes", dest="average_genes_file", type=str, help="Average of genes across cells (after scaling)")
    parser.add_argument("--gene_z", dest="zscores_gene_file", type=str,  help="z-scores for each gene in each cell")
    parser.add_argument("--probe_z", dest="zscores_probe_file", type=str,  help="z-scores for each probe in each cell")
    parser.add_argument("--corr", dest="corr_file", type=str,  help="Correlation plot")
    parser.add_argument("--problematic_probes", dest="problematic_probes_file", type=str, help="Problematic probes file")
    # parser.add_argument("--collapsed", dest="collapsed", type=mip_tools.argparse_boolean,  help="Whether the files are collapsed or not (for reporting stats)")
    parser.add_argument("--collapsed", dest="collapsed", type=str,  help="Whether the files are collapsed or not (for reporting stats)")
    parser.add_argument("--bc_stats", dest="barcodes_stats_file", type=mip_tools.argparse_bc_stats, help="path to barcode stats file, or 'None' if uncollapsed")
    parser.add_argument("--excelcounts", dest="excel_counts_file", type=str,  help="Counts per gene per cell in excel format for Jamie")
    parser.add_argument("--excelprobes", dest="excel_probes_file", type=str,  help="Summed counts per probe in excel format for Jamie")
    parser.add_argument("--excelgenes", dest="excel_gene_counts_file", type=str,  help="Summed counts per gene in excel format for Jamie")
    parser.add_argument("--excelzscores", dest="excel_zscores_file", type=str,  help="z-scores per cell in excel format for Jamie")
    
    args = parser.parse_args()

    with PdfPages(args.plots_file) as plots, open(args.stats_file, 'w') as stats, open(args.log_file, 'w') as log:
        analyze_count(args.flat_file, plots, log, stats,
                      args.collapsed, args.barcodes_stats_file, 
                      args.cell_x_gene_file, args.cell_x_gene_rescaled_file, args.cell_x_probe_file, args.cell_x_probe_rescaled_file,
                      args.sum_probes_file, args.sum_genes_file, args.average_probes_file, args.average_genes_file,
                      args.zscores_gene_file, args.zscores_probe_file, args.corr_file, args.problematic_probes_file,
                      args.excel_counts_file, args.excel_probes_file, args.excel_gene_counts_file, args.excel_zscores_file)
