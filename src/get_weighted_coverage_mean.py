#!/usr/bin/env python
import argparse
from itertools import groupby
import pandas as pd

# --------------------*********---------------------
# Description of the script
parser = argparse.ArgumentParser(description = "Script to get the weighted \
coverage mean of each contig in a set of Bins")
parser.add_argument('coverage', help="coverage table, output of jgi_summarize_contigs")
parser.add_argument('Bins', help="Bins names that are in the same folder")
parser.add_argument('out', help="coverage output")
args = parser.parse_args()

# --------------------*********---------------------
# read coverage file
coverage = pd.read_csv(args.coverage, sep="\t", index_col=0)

# --------------------*********---------------------
# open the Bins files and extract only the name and the length of each contig
weighted_coverage = coverage.loc["k99_46"] * 0
#weighted_coverage_final = coverage.loc["k99_46"] * 0
number = 0
bines = []
with open(args.Bins, 'r') as file:
    for bin in file:
        bines.append(bin.split("\n")[0])

for bin in bines:
    print bin
    with open(bin) as fasta:
        # parse each sequence by header: groupby(data, key)
        faiter = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))
        contigs = 0
        weighted_coverage = coverage.loc["k99_46"] * 0
        for record in faiter:
            # join sequence lines to one.
            for line in record:
                if line[0] == ">":
                    header1 = line.replace(">","")
                    header = header1.replace('\n','')
            contigs += 1
            seq = "".join(s.strip() for s in faiter.next())
            weighted_coverage = weighted_coverage + coverage.loc[header] * len(seq)
    wc = weighted_coverage / contigs
    Bin = str(bin.split('.fa')[0])
    wc = wc.rename(Bin)
    wc2 = pd.DataFrame(wc)
    if number == 0:
        weighted_coverage_final = wc2
    else:
        weighted_coverage_final = weighted_coverage_final.join(wc2)
    number += 1

weighted_coverage_final.T.to_csv(args.out, sep="\t", index_label="Bin")
