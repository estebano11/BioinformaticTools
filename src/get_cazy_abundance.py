#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
# --------------------*********---------------------
# Script's description
parser = argparse.ArgumentParser(description = "Get the ID abundance")
parser.add_argument('table', help="From hmmm table, outfmt with this order: \
hmm_id len query ...etc - I get the first and the third column in order to obtain \
each gene with the ID or name of another table")
parser.add_argument('cov', help="Coverage file with header")
parser.add_argument('out', help="Name/ID - coverage")
args = parser.parse_args()
# --------------------*********---------------------
# Read the Table
contigs = {}
bins = []
caz = []
with open(args.table, 'r') as file:
    for line in file:
        value = (line.split("\t")[0])[:2]
        key = (line.split("\n")[0]).split("\t")[1].split("_k99")[0]
        if key not in contigs:
            contigs[key] = [value]
            bins.append(key)
        else:
            contigs[key].append(value)
        if value not in caz:
            caz.append(value)
#df = pd.DataFrame(index=bins, columns=caz)
#df = df.fillna(0)
#for bin in contigs:
#    for gen in contigs[bin]:
#        df.loc[bin,gen] += 1
#df.to_csv(args.out, sep='\t', mode='a')
# --------------------*********---------------------
# Read coverage table
coverage = pd.read_csv(args.cov, sep="\t", index_col=0)
# # --------------------*********---------------------
# # Get for every key in table2 the total abundance
dictionary = {}
for key in contigs:
    lista = []
    list_of_values = coverage.loc["Bin.100"] * 0
    for bin in contigs[key]:
        list_of_values = list_of_values +  coverage.loc[bin]
    for i in list_of_values:
        lista.append(i)
    dictionary[key] = lista
# --------------------*********---------------------
# Print the output
with open(args.out, 'w') as handle:
    handle.write(str(dictionary))
