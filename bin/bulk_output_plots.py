#!/usr/bin/env python

## USAGE python bulk_output_plots.py <group_file> <group_plots> <group_log> <stats_file> <n_cells>
# takes in the output from umi-tools group command and does things with #s of reads, UMIs, cells, hamming distances, etc

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

def bulk_output_plots(output_plot, log, stats, groups_file, bowtie_file):
	log.write('Starting bulk analysis\n')

	# read groups file into memory
	groups = pd.read_table(groups_file, usecols=['read_id', 'umi', 'final_umi', 'unique_id'])

	# grab preliminary stats
	n_reads = len(groups)
	n_umis = groups.unique_id.max()+1

	# create column for hamming distance between the sequenced and consensus umis
	groups['hamming'] = groups.apply(lambda row: distance.hamming(row['umi'], row['final_umi']), axis=1)

	# create structure for storing how many reads support each umi cluster
	# reads_per_umi = []

	# for i, g in groups[['final_umi_count', 'unique_id']].groupby('unique_id'):
	#     reads_per_umi.append(len(g))

	reads_per_umi = groups.groupby('unique_id')['read_id'].count().values # WILL THIS WORK?

	stats.write('Reads\t{}\n'.format(n_reads))
	stats.write('UMIs\t{}\n'.format(n_umis))
	stats.write('Avg_DupRate\t{}\n'.format(n_reads / float(n_umis)))
	stats.write('Med_DupRate\t{}\n'.format(np.median(reads_per_umi)))
	stats.write('Mapping_percent\t{}\n'.format(mip_tools.grab_mapping_freq_from_bowtie(bowtie_file)))

	mip_tools.plot_reads_per_umi_bulk(reads_per_umi, output_plot)
	mip_tools.plot_umi_distance(groups.hamming, output_plot)

	log.write('Done!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Return stats from the grouped output of a hypr bulk experiment')

    parser.add_argument("-g", "--groups", dest="groups_file", type=str,  help="Groups file (from umi_tools group")
    parser.add_argument("-p", "--plots", dest="plots_file", type=str,  help="Where to output plots")
    parser.add_argument("-l", "--log", dest="log_file", type=str,  help="Where to output log")
    parser.add_argument("-s", "--stats", dest="stats_file", type=str,  help="Where to output stats")
    parser.add_argument("-b", "--bowtie", dest="bowtie_file", type=str,  help="Path to bowtie log (to grab mapping percentage from)")
    
    args = parser.parse_args()

    with PdfPages(args.plots_file) as output_plot, open(args.log_file, 'w') as log, open(args.stats_file, 'w') as stats:
    	bulk_output_plots(output_plot, log, stats, args.groups_file, args.bowtie_file)

	