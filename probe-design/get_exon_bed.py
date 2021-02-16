#!/usr/bin/env python

# technically, this doesn't actually output a bed file, but rahter a subset of the gtf (which bedtools can handle)
# turns out it was easier to do this in python than with command line tools, since there can be many matching genes
# we could try to figure out a better way of searching through the gff file..
# also, it turns out that the gtfs I have for human and mouse are formatted differently, so we need to handle them separately
# here, we simply set the gene_name_field entry for where we want to search for the gene name, and we also have check_basic to see if we need to limit our search there
# (I think that's just for mouse)
# also we have hard-coded to take the lowest transcript_id, which I think is a conservative choice but may not always be suitable...

import sys
import argparse
import pandas as pd
import openpyxl
import requests
from os import path

gtf_file = sys.argv[1]
gene_name = sys.argv[2]
out_file = sys.argv[3]
gene_name_field = sys.argv[4] # 'gene_id' # 'gene_name'
# gene_name_field = 'gene_id' # 'gene_name'
check_basic = sys.argv[5].lower() == 'true' # False for RefSeqCurated, True for CellRanger mm10 gtf (this is clunky, I know)

with open(gtf_file) as gtf:
    possible_lines = []
    for line in gtf:
        line_split = line.strip().replace('; ', '\t').split('\t')

        # determine whether to include the line
        to_include = '{} "{}"'.format(gene_name_field, gene_name) in line_split and line_split[2] == 'exon'

        # check if we have to query basic tag
        if check_basic:
            to_include = to_include and ('tag "basic"' in line_split or 'tag "basic";' in line_split)
            
        # build dict from line
        if to_include:
            temp_dict = {}
            temp_dict['chr'] = line_split[0]
            temp_dict['version'] = line_split[1]
            temp_dict['type'] = line_split[2]
            temp_dict['start'] = int(line_split[3])
            temp_dict['end'] = int(line_split[4])
            temp_dict['.'] = line_split[5]
            temp_dict['strand'] = line_split[6]
            temp_dict['.2'] = line_split[7]
            for i in range(8, len(line_split)):
                temp_dict[line_split[i].split(' ')[0]] = line_split[i].split(' ')[1]
            possible_lines.append(temp_dict) 

    # assemble dict from lines
    df = pd.DataFrame(possible_lines)

    # choose isoform/assembly with lowest transcript id
    if len(df['transcript_id'].unique()) > 1:
        tid = sorted(df['transcript_id'].unique())[0]
        df = df.loc[df.transcript_id == tid]

# reverse the order of the df if we are dealing with antisense transcript
# also looks like different gtfs are sorted in different ways
# to enforce the correct behavior, we will simply check what sense the gene is and sort ascending or descending by start, accordingly
sense_strand = df.iloc[0]['strand'] == '+'
df = df.sort_values('start', ascending=sense_strand)

# write to outfile
df[['chr', 'version', 'type', 'start', 'end', '.', 'strand', '.2', gene_name_field]].to_csv(out_file, sep='\t', header=False, index=False)

