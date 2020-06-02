#!/usr/bin/env python

import sys
import pandas as pd
import openpyxl
import math
from collections import Counter

probe_csv = sys.argv[1]
probes_for_idt = sys.argv[2]
probe_stats = sys.argv[3]

# open output
df = pd.read_csv(probe_csv, header=None)

# reformat
n_genes = len(df.columns)
n_probes = len(df) // 4
output = pd.DataFrame(columns=['name', 'seq'], index=range(2*n_probes*n_genes))

probes_per_gene = []

probe_num = 0

def get_name_from_probe(probe):
    return '-'.join(probe.split('-')[:-2])

for c in df.columns:
    for i in range(n_probes):
        probe_name = df.loc[4*i,c]
        if not pd.isnull(probe_name):
            probes_per_gene.append(get_name_from_probe(probe_name))
            probe_seq = df.loc[4*i+1,c]
            three_prime = df.loc[4*i+2,c]
            five_prime = df.loc[4*i+3,c]
            output.loc[probe_num, 'name'] = probe_name + '_5prime'
            output.loc[probe_num, 'seq'] = five_prime
            output.loc[probe_num+n_genes*n_probes, 'name'] = probe_name + '_3prime'
            output.loc[probe_num+n_genes*n_probes, 'seq'] = three_prime
            probe_num += 1

probes_per_gene = Counter(probes_per_gene)

# write to output
output.dropna().to_csv(probes_for_idt, header=False, index=False, sep='\t')

# write gene stats
with open(probe_stats, 'w') as stats:
    for gene, count in probes_per_gene.most_common():
        stats.write('{}\t{}\n'.format(gene, count))

