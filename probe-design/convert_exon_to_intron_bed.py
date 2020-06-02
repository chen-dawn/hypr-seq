#!/usr/bin/env python

import sys
import argparse
import pandas as pd
import openpyxl
import requests
from os import path

exon_bed = sys.argv[1]
intron_bed = sys.argv[2]

distance_from_tss = 5000

# read exon bed
exons = pd.read_table(exon_bed, header=None)

# intialize intron df (copy of all the columns, minus the last row)
introns = exons.copy().loc[range(len(exons)-1)]

# replace values with intron coords
for r in range(len(exons)-1):
    introns.loc[r, [2, 3, 4]] = ['intron', exons.loc[r, 4], exons.loc[r+1,3]]

introns[3] = introns[3] + 1
introns[4] = introns[4] - 1

# figure out the strandedness
if '-' in introns[6].values: # negative
    # get first introns first
    introns = introns.sort_values(3, ascending=False)

    # get TSS
    tss = exons.loc[len(exons)-1, 4]

    # get only things within some distance (5kb) of tss
    introns_to_use = introns.loc[tss - introns[4] < distance_from_tss].copy()
    introns_to_use.index = range(len(introns_to_use))
    introns_to_use.loc[len(introns_to_use)-1, 3] = max(introns_to_use.loc[len(introns_to_use)-1, 3], tss-5000)
else: # positive
    # get TSS
    tss = exons.loc[0, 3]

    # get only things within some distance (5kb) of tss
    introns_to_use = introns.loc[introns[3] - tss < distance_from_tss].copy()
    introns_to_use.loc[len(introns_to_use)-1, 4] = min(introns_to_use.loc[len(introns_to_use)-1, 4], tss+5000)

introns_to_use.to_csv(intron_bed, sep='\t', header=False, index=False)
