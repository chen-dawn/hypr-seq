#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
import argparse

import mip_tools

def collate(output_file, input_files):
    with open(output_file, 'w') as output:
        output.write('Done!')
        for file in input_files:
            with open(file, 'r') as f:
                pass
                # output.write('{}\n'.format(len(f.readlines())))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tmp file to collate all the outputs and make sure we have them.')
    
    parser.add_argument("-o", "--output", dest="output_file", type=str, help="Where to store the output table")
    parser.add_argument("-i", "--inputs", dest="input_files", type=str, nargs='+', help="Input count files")
    
    args = parser.parse_args()

    collate(args.output_file, args.input_files)
