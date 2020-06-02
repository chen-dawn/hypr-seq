#!/usr/bin/env python

import sys
import argparse
import pandas as pd
import openpyxl
import requests
from os import path

gtf_file = sys.argv[1]
gene_name = sys.argv[2]
out_file = sys.argv[3]

with open(gtf_file) as gtf:
    possible_lines = []
    for line in gtf:
        line_split = line.strip().replace('; ', '\t').split('\t')
        if 'gene_name "{}"'.format(gene_name) in line_split and line_split[2] == 'exon' and 'tag "basic"' in line_split:
            temp_dict = {}
            temp_dict['chr'] = line_split[0]
            temp_dict['version'] = line_split[1]
            temp_dict['type'] = line_split[2]
            temp_dict['start'] = line_split[3]
            temp_dict['end'] = line_split[4]
            temp_dict['.'] = line_split[5]
            temp_dict['strand'] = line_split[6]
            temp_dict['.2'] = line_split[7]
            for i in range(8, len(line_split)):
                temp_dict[line_split[i].split(' ')[0]] = line_split[i].split(' ')[1]
            possible_lines.append(temp_dict) 
    df = pd.DataFrame(possible_lines)
    if len(df['transcript_id'].unique()) > 1:
        tid = sorted(df['transcript_id'].unique())[0]
        df = df.loc[df.transcript_id == tid]

df[['chr', 'version', 'type', 'start', 'end', '.', 'strand', '.2', 'gene_name', 'exon_number']].sort_values('exon_number').to_csv(out_file, sep='\t', header=False, index=False)
    




