#!/usr/bin/env python

# usage: ./import_mips_new.py infile outfile

import pandas as pd
import argparse

def read_probes_to_fa(in_file, out_file):
    # check file format, excel or txt
    if '.xls' in in_file:
        # excel file
        file = pd.read_excel(in_file, header=None, names=['name', 'seq'])
    else:
        # txt file
        file = pd.read_table(in_file, header=None, names=['name', 'seq'])
    
    for probe in file.itertuples():
        probe_name, sequence = probe.name, probe.seq
        out_file.write('>{}\n'.format(probe_name.replace(' ', '')))
        out_file.write(sequence)
        out_file.write('\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert probes in excel format to a fasta file for bowtie-build')
    
    parser.add_argument("-i", dest="in_file", type=str,  help="Input probe file (excel or txt")
    parser.add_argument("-o", dest="out_file", type=str,  help="Output fasta file")

    args = parser.parse_args()

    with open(args.out_file, 'w') as out_file:
        read_probes_to_fa(args.in_file, out_file)
