#!/usr/bin/env python

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
from scipy.stats import ttest_ind
import os
from os import path

# ncs = ['Brainbar001', 'Brainbar003', 'SplintRbar002', 'SplintRbar23', 'SplintRbar26', 'SplintRbar27', 'SplintRbar28', 'SplintRbar30']
# bc_to_guide = {'Brainbar001': 'NC-1', 'Brainbar002': 'e-GATA1-1', 'Brainbar003': 'NC-2', 'Brainbar004': 'e-GATA1-2', 'Brainbar005': 'e-HDAC6-1', 'Brainbar007': 'e-HDAC6-2', 'Brainbar008': 'HDAC6-TSS-1', 'SplintRbar001': 'GATA1-TSS-1', 'SplintRbar002': 'NC-5', 'SplintRbar4': 'GATA1-TSS-2', 'SplintRbar23': 'NC-4', 'SplintRbar26': 'NC-6', 'SplintRbar27': 'NC-7', 'SplintRbar28': 'NC-8', 'SplintRbar30': 'NC-9', 'SplintRbar5': 'HDAC6-TSS-2'}
# ls = ['NC-1','NC-2','NC-4','NC-5','NC-6','NC-7','NC-8','NC-9','GATA1-TSS-1','GATA1-TSS-2','HDAC6-TSS-1','HDAC6-TSS-2','e-GATA1-1','e-GATA1-2','e-HDAC6-1','e-HDAC6-2']

def analyze_p16(collapsed_counts_file, uncollapsed_counts_file, barcode_stats_file, 
                guide_info_file, target_genes, housekeeping_genes, 
                plots, log, stats, barcodes_file, chastity, entropy_thresh, 
                guide_stats_file, element_stats_file, filter_dcas9, dcas9_transcript):

    # read in counts
    collapsed_counts = pd.read_table(collapsed_counts_file, index_col=0)
    uncollapsed_counts = pd.read_table(uncollapsed_counts_file, index_col=0)

    # compute entropy
    entropy_collapsed = mip_tools.compute_entropy(collapsed_counts)
    entropy_uncollapsed = mip_tools.compute_entropy(uncollapsed_counts)

    # if we have barode classes, plot those
    if barcode_stats_file:
        # grab barcodes corresponding to singletons and multiplets
        barcode_stats = pd.read_table(barcode_stats_file, index_col=0)

        singletons = barcode_stats.loc[barcode_stats['Singleton']==True].index.tolist()
        multiplets = barcode_stats.loc[barcode_stats['Multiplet']==True].index.tolist()
        giant = barcode_stats.loc[barcode_stats['Giant']==True].index.tolist()

        # plot entropy distributions separately
        hist, bins = np.histogram(entropy_uncollapsed.dropna(), bins='auto')

        plt.hist(entropy_uncollapsed.loc[multiplets].dropna(), bins=bins, alpha=0.5, label='Multiplet')
        plt.hist(entropy_uncollapsed.loc[singletons].dropna(), bins=bins, alpha=0.5, label='Singletons')
        plt.hist(entropy_uncollapsed.loc[giant].dropna(), bins=bins, alpha=0.5, label='Giant')

        plt.legend()


    # common plotting stuff
    plt.xlabel('Entropy (bits)')
    plt.ylabel('Number of barcodes')
    plt.title('Entropy plot (all barcodes)')
    plt.tight_layout()
    plots.savefig()
    plt.close() 

    plt.hist(entropy_collapsed.loc[entropy_collapsed.index.intersection(singletons)].dropna(), bins='auto')
    plt.xlabel('Entropy (bits)')
    plt.ylabel('Number of barcodes')
    plt.title('Entropy plot (collapsed barcodes)')
    plt.tight_layout()
    plots.savefig()
    plt.close()

    collapsed_counts = collapsed_counts.loc[collapsed_counts.index.intersection(singletons+multiplets)]

    # assign barcodes to cells
    unassigned_rate_uncollapsed = mip_tools.assign_barcode(uncollapsed_counts, chastity)
    stats.write('Correctly Assigned Cells (uncollapsed)\t{:.2f}%\n'.format((1-unassigned_rate_uncollapsed)*100))
    stats.write('Unassigned Cells (uncollapsed)\t{:.2f}%\n'.format(unassigned_rate_uncollapsed*100))

    unassigned_rate_collapsed = mip_tools.assign_barcode(collapsed_counts, chastity)
    stats.write('Correctly Assigned Cells (collapsed)\t{:.2f}%\n'.format((1-unassigned_rate_collapsed)*100))
    stats.write('Unassigned Cells (collapsed)\t{:.2f}%\n'.format(unassigned_rate_collapsed*100))

    # compute the cell lambda from p16
    cell_lambda = mip_tools.calculate_cell_lambda_from_p16_counts(collapsed_counts)
    stats.write('cell_lambda\t{:.3f}\n'.format(cell_lambda))

    # compute mixing fraction
    mix_frac_uncollapsed = (entropy_uncollapsed > entropy_thresh).mean()
    stats.write('Mixing (uncollasped entropy > {:.2f})\t{:.2f}%\n'.format(entropy_thresh, mix_frac_uncollapsed*100))

    mix_frac_collapsed = (entropy_collapsed > entropy_thresh).mean()
    stats.write('Mixing (collapsed entropy > {:.2f})\t{:.2f}%\n'.format(entropy_thresh, mix_frac_collapsed*100))

    # write barcodes to file
    collapsed_counts[['Barcode']].to_csv(barcodes_file, sep='\t')

    # use barcodes to annotate the rows with additional info
    guide_info = pd.read_table(guide_info_file)
    collapsed_counts = collapsed_counts.merge(guide_info, left_on='Barcode', right_on='barcode')

    # what to do about BFP+ cells
    # if there's a BFP column, take only the positive cells
    if filter_dcas9:
        bfp_col = collapsed_counts.filter(like=dcas9_transcript).columns[0]
        min_bfp = mip_tools.get_and_plot_bfp_threshold(collapsed_counts, plots, BFP_col=bfp_col)
        stats.write('BFP threshold\t{:.2f}\n'.format(min_bfp))
        stats.write('Percent cells with BFP below threshold\t{:.2f}%\n'.format(100 * sum((collapsed_counts.Barcode != 'Unassigned') & (collapsed_counts[bfp_col] < min_bfp)) / sum(collapsed_counts.Barcode != 'Unassigned')))
        collapsed_counts = collapsed_counts.loc[collapsed_counts[bfp_col] >= min_bfp]

    mip_tools.analyze_p16_data(collapsed_counts, housekeeping_genes, target_genes, plots, guide_stats_file, element_stats_file, log)

    log.write('Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze the count matrix from a Hypr-seq experiment and return useful plots.')
    
    parser.add_argument("--collapsed", dest="collapsed", type=str,  help="Collapsed counts file")
    parser.add_argument("--uncollapsed", dest="uncollapsed", type=str,  help="Uncollapsed counts file")
    parser.add_argument("--barcode_stats", dest="barcode_stats", type=str,  help="Barcode multiplicity stats file")
    parser.add_argument("--guide_info", dest="guide_info_file", type=str,  help="Information about guide+element status")
    parser.add_argument("--target_genes", dest="target_gene_list", type=str, default='GATA1' help="Comma-separated list of genes to plot knockdown of")
    parser.add_argument("--housekeeping_genes", dest="housekeeping_gene_list", type=str, default='GAPDH' help="Comma-separated list of housekeeping genes to normalize by")
    parser.add_argument("--chastity", dest="chastity", type=float, default=6., help="How strict to be for barcode assignment")
    parser.add_argument("--entropy_thresh", dest="entropy_thresh", type=float, default=0.5, help="Entropy cutoff to call mixed cell")
    parser.add_argument("--outdir", dest="outdir", type=str, help="Output directory for results")
    parser.add_argument("--name", dest="sample_name", type=str, help="Name for outputting plots, stats, etc.")
    parser.add_argument("--codedir", dest="codedir", type=str, help="Path to hypr code directory")
    parser.add_argument("--filter_dcas9", dest="filter_dcas9", action='store_true', help="Whether to filter cells by a marker transcript for dCas9 expression")
    parser.add_argument("--dcas9_transcript", dest="dcas9_transcript", type=str, default='BFP', help="Trasncript to filter on to get high dCas9-expressing cells")

    args = parser.parse_args()

    ## IMPORT MIP TOOLS FROM DIFFERENT DIR
    sys.path.insert(0, '{}/bin'.format(args.codedir))
    import mip_tools

    # format output file names
    plots_file = path.join(args.outdir, '{}.crispr_screen.plots.pdf'.format(args.sample_name))
    log_file = path.join(args.outdir, '{}.crispr_screen.log.txt'.format(args.sample_name))
    stats_file = path.join(args.outdir, '{}.crispr_screen.stats.txt'.format(args.sample_name))
    barcodes_file = path.join(args.outdir, '{}.crispr_screen.barcode_assignments.txt'.format(args.sample_name))
    guide_stats_file = path.join(args.outdir, '{}.crispr_screen.guide_stats.txt'.format(args.sample_name))
    element_stats_file = path.join(args.outdir, '{}.crispr_screen.element_stats.txt'.format(args.sample_name))

    target_genes = args.target_gene_list.split(',')
    housekeeping_genes = args.housekeeping_gene_list.split(',')

    with PdfPages(plots_file) as plots, open(stats_file, 'w') as stats, open(log_file, 'w') as log:
        analyze_p16(args.collapsed, args.uncollapsed, args.barcode_stats, args.guide_info_file, target_genes, housekeeping_genes, plots, log, stats, barcodes_file, args.chastity, args.entropy_thresh, guide_stats_file, element_stats_file, args.filter_dcas9, args.dcas9_transcript)
        
