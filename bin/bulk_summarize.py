#!/usr/bin/env python

import sys
from os import path
import pandas as pd
import numpy as np
import argparse

import mip_tools

def summarize_bulk(output_stats_file, input_stats_files, output_counts_file, input_counts_files):
    input_counts = [pd.read_table(f, header=None, index_col=0, names=[path.basename(f).split('_bulk.counts')[0]]) for f in input_counts_files]
    merged_counts = pd.concat(input_counts, axis=1).fillna(0)
    merged_counts.to_csv(output_counts_file, sep='\t')

    input_stats = [pd.read_table(f, header=None, index_col=0, names=[path.basename(f).split('_bulk.stats')[0]]) for f in input_stats_files]
    merged_stats = pd.concat(input_stats, axis=1).fillna(0)
    merged_stats.to_csv(output_stats_file, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze the count matrix from a Hypr-seq experiment and return useful plots.')
    
    parser.add_argument("-t", "--output_stats", dest="output_stats_file", type=str, help="Where to store the stats table")
    parser.add_argument("-s", "--input_stats", dest="input_stats_files", type=str, nargs='+', help="Input stats files")
    parser.add_argument("-d", "--output_counts", dest="output_counts_file", type=str, help="Where to store the output table")
    parser.add_argument("-c", "--input_counts", dest="input_counts_files", type=str, nargs='+', help="Input count files")
    
    args = parser.parse_args()

    summarize_bulk(args.output_stats_file, args.input_stats_files, args.output_counts_file, args.input_counts_files)
