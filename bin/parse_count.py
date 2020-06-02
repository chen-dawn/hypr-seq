#!/usr/bin/env python

## USAGE python parse_count.py <raw_count_file> <collapsed_count_file> <count_plots> 
# takes the raw counts (per individual probe) and collapses to probes sharing the same transcript 

import sys
import mip_tools
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')
import pandas as pd
import numpy as np
import argparse

plot_scale = 3

def parse_count(in_counts_file, out_counts_file, out_probes_file, output_plot):
	# read input file
	in_counts = pd.read_table(in_counts_file, sep='\\s+', header=None, names=['Counts', 'LongName'])

	# process to add new columns
	in_counts['BaseName'] = in_counts.apply(lambda row: mip_tools.extract_base_probe_name(row.LongName), axis=1)
	in_counts['ProbeNumber'] = in_counts.apply(lambda row: mip_tools.extract_probe_number(row.LongName), axis=1)
	in_counts = in_counts.groupby(['BaseName', 'ProbeNumber']).agg(np.sum).reset_index()

	# print out plots per transcript
	for probe, subset in in_counts.groupby('BaseName'):
		subset.plot.bar(x='ProbeNumber', y='Counts')
		plt.title('Counts per probe: {}'.format(probe))
		plt.yscale('log', nonposy='clip')
		plt.ylim(1, plot_scale*subset.Counts.values.max())
		plt.tight_layout()
		output_plot.savefig()

	# pivot and save to file
	pivoted = in_counts.pivot(index='BaseName', columns='ProbeNumber', values='Counts')
	pivoted.to_csv(out_probes_file, sep='\t')

	# collapse to transcript counts
	collapsed = in_counts.groupby('BaseName')['Counts'].sum()

	# plot bulk transcript counts
	collapsed.plot.bar()
	plt.title('Probe Counts')
	plt.yscale('log', nonposy='clip')
	plt.ylim(1, plot_scale*collapsed.values.max())
	output_plot.savefig()

	# write output
	collapsed.to_csv(out_counts_file, sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze the count matrix from a Hypr-seq experiment and return useful plots.')
    
    parser.add_argument("-r", "--raw", dest="in_counts_file", type=str,  help="Raw counts file (flat and not filtered)")
    parser.add_argument("-o", "--plots", dest="out_plot_file", type=str,  help="Where to output plots")
    parser.add_argument("-c", "--counts", dest="out_counts_file", type=str,  help="Output counts file")
    parser.add_argument("-p", "--probes", dest="out_probes_file", type=str,  help="Output probes file")
    
    args = parser.parse_args()

    with PdfPages(args.out_plot_file) as output_plot:
    	parse_count(args.in_counts_file, args.out_counts_file, args.out_probes_file, output_plot)

		