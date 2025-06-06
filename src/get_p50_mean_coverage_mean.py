#!/usr/bin/env python3
import argparse
from itertools import groupby
import pandas as pd

# --------------------*********---------------------
# Description of the script
parser = argparse.ArgumentParser(description = "Script to get the mean \
coverage of each Bins")
parser.add_argument('coverage', help="coverage table, output of jgi_summarize_contigs")
parser.add_argument('Bins', help="Bins names that are in the same folder")
parser.add_argument('n50', help="File with Bin name and N50")
parser.add_argument('out', help="coverage output")
args = parser.parse_args()

# --------------------*********---------------------
# read N50 file
n50 = []
with open(args.n50, 'r') as file:
    for bin in file:
        n50.append([bin.split("\t")[0],int(bin.split("\t")[1].split("\n")[0])])

# read coverage file
coverage = pd.read_csv(args.coverage, sep="\t", index_col=0)

# --------------------*********---------------------
# open the Bins files and extract only the name and the length of each contig
mean_coverage = coverage.loc["Bin.100_k99_4421153"] * 0 # ver una mejor forma de hacer esto
number = 0
bines = []
with open(args.Bins, 'r') as file:
    for bin in file:
        bines.append(bin.split("\n")[0])
def Grep(bin):
    for i in n50:
        if i[0] == bin:
            return i[1]

for bin in bines:
    N50 = Grep(bin)
    print bin
    with open(bin) as fasta:
        # parse each sequence by header: groupby(data, key)
        faiter = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))
        contigs = 0
        size = 0
        mean_coverage = coverage.loc["Bin.100_k99_4421153"] * 0
        for record in faiter:
            # join sequence lines to one.
            for line in record:
                if line[0] == ">":
                    header1 = line.replace(">","")
                    header = header1.replace('\n','')
            seq = "".join(s.strip() for s in faiter.next())
            if len(seq) > N50:
                mean_coverage = mean_coverage + coverage.loc[header]
                contigs += 1
                size += len(seq)
    wc = mean_coverage / contigs
    Bin = str(bin.split('.fa')[0])
    print (wc)
    wc = wc.rename(Bin)
    wc2 = pd.DataFrame(wc)
    if number == 0:
        mean_coverage_final = wc2
    else:
        mean_coverage_final = mean_coverage_final.join(wc2)
    number += 1

mean_coverage_final.T.to_csv(args.out, sep="\t", index_label="Bin")
