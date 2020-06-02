#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')
import itertools
import distance
import argparse

import mip_tools

def summarize_stats_from_groups_file(groups):
    reads_per_umi = groups.groupby('unique_id')['read_id'].count().values
    reads_per_cell = groups.groupby('CBC')['unique_id'].count().values
    umis_per_cell = groups.groupby('CBC')['unique_id'].agg(lambda df: len(df.unique())).values
    total_reads = len(groups)
    total_umis = len(groups['unique_id'].unique())
    total_cells = len(groups['CBC'].unique())
    return {'reads/umi': reads_per_umi,
            'reads/cell': reads_per_cell,
            'umis/cell': umis_per_cell, 
            'reads': total_reads, 
            'umis': total_umis, 
            'cells': total_cells}

def write_stats_to_file(group_stats, stats, stats_suffix, bowtie_file):
    stats.write('Reads{}\t{}\n'.format(stats_suffix, group_stats['reads']))
    stats.write('UMIs{}\t{}\n'.format(stats_suffix, group_stats['umis']))
    stats.write('Cells{}\t{}\n'.format(stats_suffix, group_stats['cells']))

    stats.write('Avg_DupRate{}\t{}\n'.format(stats_suffix, np.mean(group_stats['reads/umi'])))
    stats.write('Med_DupRate{}\t{}\n'.format(stats_suffix, np.median(group_stats['reads/umi'])))
    stats.write('Sequencing_saturation{}\t{}\n'.format(stats_suffix, 1 - group_stats['umis'] / group_stats['reads']))

    stats.write('Avg_UMIs/cell{}\t{}\n'.format(stats_suffix, np.mean(group_stats['umis/cell'])))
    stats.write('Med_UMIs/cell{}\t{}\n'.format(stats_suffix, np.median(group_stats['umis/cell'])))

    stats.write('Avg_Reads/cell{}\t{}\n'.format(stats_suffix, np.mean(group_stats['reads/cell'])))
    stats.write('Med_Reads/cell{}\t{}\n'.format(stats_suffix, np.median(group_stats['reads/cell'])))

    stats.write('Mapping_percent{}\t{}\n'.format(stats_suffix, mip_tools.grab_mapping_freq_from_bowtie(bowtie_file)))

def plot_downsample(fracs, total_umis, reads_per_umi, saturations, output_plot):
    plt.plot(fracs, total_umis)
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlabel('Fraction of reads')
    plt.ylabel('Number of UMIs')
    plt.title('Downsampled UMIs')
    plt.tight_layout()
    output_plot.savefig()
    plt.close()

    plt.plot(fracs, reads_per_umi)
    # plt.xscale('log')
    plt.xlabel('Fraction of reads')
    plt.ylabel('Duplication rate')
    plt.title('Downsampled duplication rate')
    plt.tight_layout()
    output_plot.savefig()
    plt.close()

    plt.plot(fracs, saturations)
    # plt.xscale('log')
    plt.xlabel('Fraction of reads')
    plt.ylabel('Sequencing saturation')
    plt.title('Downsampled saturation rate')
    plt.tight_layout()
    output_plot.savefig()
    plt.close()

dup_fracs = np.linspace(1e-5, 1., num=25) # [1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001]

def downsample_reads(groups, output_plot, fracs=dup_fracs):
    total_umis, reads_per_umi, saturations = [], [], []

    for f in fracs:
        gs = summarize_stats_from_groups_file(groups.sample(frac=f))
        total_umis.append(gs['umis'])
        reads_per_umi.append(gs['reads/umi'].mean())
        saturations.append(1 - gs['umis'] / gs['reads'])

    plot_downsample(fracs, total_umis, reads_per_umi, saturations, output_plot)


def analyze_emul_duplication(output_plot, log, stats, overlap_thresh, groups_file, bowtie_file, collapsed):
    log.write('Starting emul analysis\n')
    # stats_suffix = '_collapsed' if collapsed else '_uncollapsed'
    stats_suffix = '_'+collapsed

    # read groups file into memory
    groups = pd.read_table(groups_file, usecols=['read_id', 'umi', 'final_umi', 'unique_id'])    

    # add column for cell barcode info 
    groups['CBC'] = groups.apply(lambda row: mip_tools.extract_cbc(row['read_id']), axis=1)

    # write stats
    group_stats = summarize_stats_from_groups_file(groups)
    write_stats_to_file(group_stats, stats, stats_suffix, bowtie_file)

    # make output plots
    mip_tools.plot_histogram_and_cdf(group_stats['reads/umi'], output_plot, 'Reads per UMI', 'UMIs', 'Reads per UMI')
    mip_tools.plot_histogram_and_cdf(group_stats['reads/cell'], output_plot, 'Reads per Cell', 'Cells', 'Reads per Cell')
    mip_tools.plot_histogram_and_cdf(group_stats['umis/cell'], output_plot, 'UMIs per Cell', 'Cells', 'UMIs per Cell')

    # # create column for hamming distance between the sequenced and consensus umis
    # groups['hamming'] = groups.apply(lambda row: distance.hamming(row['umi'], row['final_umi']), axis=1)
    # mip_tools.plot_umi_distance(groups.hamming, output_plot)

    # downsampling
    downsample_reads(groups, output_plot)


    log.write('Done!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Return stats from the grouped output of a hypr emulsion experiment')

    parser.add_argument("-g", "--groups", dest="groups_file", type=str,  help="Groups file (from umi_tools group")
    parser.add_argument("-p", "--plots", dest="plots_file", type=str,  help="Where to output plots")
    parser.add_argument("-l", "--log", dest="log_file", type=str,  help="Where to output log")
    parser.add_argument("-s", "--stats", dest="stats_file", type=str,  help="Where to output stats")
    parser.add_argument("-b", "--bowtie", dest="bowtie_file", type=str,  help="Path to bowtie log (to grab mapping percentage from)")
    parser.add_argument("-t", "--threshold", dest="threshold", type=int,  help="Number of shared UMIs to determine beads came from same cell", default=10)
    # parser.add_argument("-u", "--collapsed", dest="collapsed", type=mip_tools.argparse_boolean,  help="Whether the files are collapsed or not (for reporting stats)")
    parser.add_argument("-u", "--collapsed", dest="collapsed", type=str,  help="Whether the files are collapsed or not (for reporting stats)")


    args = parser.parse_args()

    with PdfPages(args.plots_file) as output_plot, open(args.log_file, 'w') as log, open(args.stats_file, 'w') as stats:
        analyze_emul_duplication(output_plot, log, stats, args.threshold, args.groups_file, args.bowtie_file, args.collapsed)
