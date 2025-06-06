#!/usr/bin/env python3
import argparse
import pandas as pd
import csv
# --------------------*********---------------------
# Script's description
parser = argparse.ArgumentParser(description = "From bin-class table and caz-Bin  \
and coverage table, get the caz-abundance by class")
parser.add_argument('table', help="From hmmm table, outfmt with this order: \
hmm_id query - I get the first and the third column in order to obtain \
each gene with the ID or name of another table")
parser.add_argument('clase', help="Class table: Bin class tab separated with no header")
parser.add_argument('cov', help="Coverage file with header")
parser.add_argument('out', help="caz-coverage")
args = parser.parse_args()
# --------------------*********---------------------
# Read the Table
contigs = {}
with open(args.table, 'r') as file:
    for line in file:
        #key = line.split("\t")[0]
        key = (line.split("\t")[0])
        #value = (line.split("\t")[2]).split("_k99")[0]
        value = (line.split("\n")[0]).split("\t")[1].split("_k99")[0]
        if key not in contigs:
            contigs[key] = [value]
        else:
            contigs[key].append(value)
clase = {}
max = 0
with open(args.clase, 'r') as file:
    for line in file:
        key = line.split("\t")[0]
        value = int((line.split("\n")[0]).split("\t")[1])
        clase[key] = value
        if value > max:
            max = value
# --------------------*********---------------------
# Read coverage table
coverage = pd.read_csv(args.cov, sep="\t", index_col=0)
#coverage[coverage>=1] = 1
#coverage[coverage<1] = 0
#coverage.to_csv("depth.tmp", sep='\t', mode='a')
# --------------------*********---------------------
m = 0
while m <= max:
    if m == 1 or m == 5:
        modulo = "modulo_1-5"
        subset = {}
        for key in contigs:
            for bin in contigs[key]:
                if bin in clase and clase[bin] == m:
                    if key not in subset:
                        subset[key] = [bin]
                    else:
                        subset[key].append(bin)
        dictionary = {}
        for key in subset:
            l = []
            list_of_values = coverage.loc["Bin.100"] * 0
            for bin in subset[key]:
                list_of_values = list_of_values +  coverage.loc[bin]
            for i in list_of_values:
                l.append(i)
            dictionary[key] = l
        w = csv.writer(open(args.out + "_" + modulo + ".csv", 'w'))
        for key, val in dictionary.items():
                w.writerow([key,val])
        m += 1
