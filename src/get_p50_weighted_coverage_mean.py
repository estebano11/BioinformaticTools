#!/usr/mag/env python
import argparse
from itertools import groupby
import pandas as pd

# --------------------*********---------------------
# Description of the script
parser = argparse.ArgumentParser(description = "Script to get the weighted \
coverage mean of each contig in a set of Bins")
parser.add_argument('coverage', help="coverage table, output of jgi_summarize_contigs")
parser.add_argument('Bins', help="Bins names that are in the same folder")
parser.add_argument('n50', help="File with Bin name and N50")
parser.add_argument('out', help="coverage output")
args = parser.parse_args()

# --------------------*********---------------------
# read N50 file
n50 = []
with open(args.n50, 'r') as file:
    for mag in file:
        n50.append([mag.split("\t")[0],int(mag.split("\t")[1].split("\n")[0])])
# read coverage file
coverage = pd.read_csv(args.coverage, sep="\t", index_col=0)
lista_contigs_cov = list(coverage.index)
# --------------------*********---------------------
# open the Bins files and extract only the name and the length of each contig
# Extraigo el primer elemento de la lista de cobertura para armar mi dataframe
weighted_coverage = coverage.loc[lista_contigs_cov[0]] * 0
number = 0
mags = []
with open(args.Bins, 'r') as file:
    for mag in file:
        mags.append(mag.split("\n")[0])
def Grep(mag):
    for i in n50:
        if i[0] == mag:
            return i[1]
for mag in mags:
    print (mag)
    N50 = Grep(mag)
    print (type(N50), N50)
    with open(mag) as fasta:
        # parse each sequence by header: groupby(data, key)
        faiter = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))
        contigs = 0
        size = 0
        weighted_coverage = coverage.loc[lista_contigs_cov[0]] * 0
        for record in faiter:
            # join sequence lines to one.
            for line in record:
                if line[0] == ">":
                    header1 = line.replace(">","")
                    header = header1.replace('\n','')
            seq = "".join(s.strip() for s in faiter.__next__())
            if len(seq) > N50:
                weighted_coverage = weighted_coverage + coverage.loc[header] * len(seq)
                contigs += 1
                size += len(seq)
    #wc = weighted_coverage / (contigs * size) # No debería multiplciar por los contigs, ya que con
    #                                           tamaño del genoma estaría
    wc = weighted_coverage / (size)
    Bin = str(mag.split('.fa')[0])
    wc = wc.rename(Bin)
    wc2 = pd.DataFrame(wc)
    if number == 0:
        weighted_coverage_final = wc2
    else:
        weighted_coverage_final = weighted_coverage_final.join(wc2)
    number += 1

weighted_coverage_final.T.to_csv(args.out, sep="\t", index_label="Bin")
