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
from scipy.stats import ttest_ind

import mip_tools

ncs = ['Brainbar001', 'Brainbar003', 'SplintRbar002', 'SplintRbar23', 'SplintRbar26', 'SplintRbar27', 'SplintRbar28', 'SplintRbar30']
bc_to_guide = {'Brainbar001': 'NC-1', 'Brainbar002': 'e-GATA1-1', 'Brainbar003': 'NC-2', 'Brainbar004': 'e-GATA1-2', 'Brainbar005': 'e-HDAC6-1', 'Brainbar007': 'e-HDAC6-2', 'Brainbar008': 'HDAC6-TSS-1', 'SplintRbar001': 'GATA1-TSS-1', 'SplintRbar002': 'NC-5', 'SplintRbar4': 'GATA1-TSS-2', 'SplintRbar23': 'NC-4', 'SplintRbar26': 'NC-6', 'SplintRbar27': 'NC-7', 'SplintRbar28': 'NC-8', 'SplintRbar30': 'NC-9', 'SplintRbar5': 'HDAC6-TSS-2'}
ls = ['NC-1','NC-2','NC-4','NC-5','NC-6','NC-7','NC-8','NC-9','GATA1-TSS-1','GATA1-TSS-2','HDAC6-TSS-1','HDAC6-TSS-2','e-GATA1-1','e-GATA1-2','e-HDAC6-1','e-HDAC6-2']

def analyze_p16(collapsed_counts_file, uncollapsed_counts_file, barcode_stats_file, plots, log, stats, barcodes_file, chastity, entropy_thresh, guide_stats_file, element_stats_file):
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

        # multiplets = pd.read_table(multiplets_file, squeeze=True)
        # singletons = pd.read_table(singletons_file, squeeze=True)        

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

    # plot two entropy distributions separately
    # plt.hist(entropy_collapsed.dropna(), bins='auto')
    # plt.hist(entropy_collapsed.loc[entropy_collapsed.index.intersection(singletons+multiplets)].dropna(), bins='auto')
    plt.hist(entropy_collapsed.loc[entropy_collapsed.index.intersection(singletons)].dropna(), bins='auto')
    plt.xlabel('Entropy (bits)')
    plt.ylabel('Number of barcodes')
    plt.title('Entropy plot (collapsed barcodes)')
    plt.tight_layout()
    plots.savefig()
    plt.close()

    # note, changed this from singletons to singletons and multiplets 191110
    # potentially in some of the noisy KD plots it was in experiments with lots of doublets that we were ignoring?
    # think this is still the optimal strategy, anyway..
    collapsed_counts = collapsed_counts.loc[collapsed_counts.index.intersection(singletons+multiplets)]
    # collapsed_counts = collapsed_counts.loc[collapsed_counts.index.intersection(singletons)]

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
    collapsed_counts['guide'] = collapsed_counts.apply(lambda row: bc_to_guide.get(row.Barcode, row.Barcode) if row.Barcode != 'Unassigned' else row.Barcode, axis=1)
    collapsed_counts['type'] = collapsed_counts.apply(lambda row: '-'.join(row['guide'].split('-')[:-1]), axis=1)
    collapsed_counts['target'] = collapsed_counts.apply(lambda row: 'negative_control' if row['type'] == 'NC' else row['type'], axis=1)

    # what to do about BFP+ cells
    # if there's a BFP column, take only the positive cells
    if not collapsed_counts.filter(like='BFP').columns.empty:
        bfp_col = collapsed_counts.filter(like='BFP').columns[0]
        min_bfp = mip_tools.get_and_plot_bfp_threshold(collapsed_counts, plots, BFP_col=bfp_col)
        stats.write('BFP threshold\t{:.2f}\n'.format(min_bfp))
        stats.write('Percent cells with BFP below threshold\t{:.2f}%\n'.format(100 * sum((collapsed_counts.Barcode != 'Unassigned') & (collapsed_counts[bfp_col] < min_bfp)) / sum(collapsed_counts.Barcode != 'Unassigned')))
        collapsed_counts = collapsed_counts.loc[collapsed_counts[bfp_col] >= min_bfp]

    mip_tools.analyze_p16_data(collapsed_counts, plots, guide_stats_file, element_stats_file, log)

    log.write('Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze the count matrix from a Hypr-seq experiment and return useful plots.')
    
    parser.add_argument("-c", "--collapsed", dest="collapsed", type=str,  help="Collapsed counts file ")
    parser.add_argument("-d", "--uncollapsed", dest="uncollapsed", type=str,  help="Uncollapsed counts file ")
    # parser.add_argument("-m", "--multiplets", dest="multiplets_file", type=str,  help="Multiplet barcodes file ")
    parser.add_argument("-r", "--barcode_stats", dest="barcode_stats", type=str,  help="Singleton barcodes file ")
    parser.add_argument("-p", "--plots", dest="plots_file", type=str,  help="Where to output plots")
    parser.add_argument("-l", "--log", dest="log_file", type=str,  help="Where to output log")
    parser.add_argument("-s", "--stats", dest="stats_file", type=str,  help="Where to output stats")
    parser.add_argument("-b", "--barcodes", dest="barcodes_file", type=str,  help="Output barcode assignments per cell")
    parser.add_argument("-t", "--chastity", dest="chastity", type=float, default=6., help="How strict to be for barcode assignment")
    parser.add_argument("-e", "--entropy_thresh", dest="entropy_thresh", type=float, default=0.5, help="Entropy cutoff to call mixed cell")
    parser.add_argument("-k", "--guides", dest="guide_stats_file", type=str,  help="Where to output knockdown info for guides")
    parser.add_argument("-v", "--elements", dest="element_stats_file", type=str,  help="Where to output knockdown info for elements")


    args = parser.parse_args()

    with PdfPages(args.plots_file) as plots, open(args.stats_file, 'w') as stats, open(args.log_file, 'w') as log:
        analyze_p16(args.collapsed, args.uncollapsed, args.barcode_stats, plots, log, stats, args.barcodes_file, args.chastity, args.entropy_thresh, args.guide_stats_file, args.element_stats_file)






 


