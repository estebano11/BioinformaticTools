#!/usr/bin/env python3
import argparse
import pandas as pd
# --------------------*********---------------------
# Script's description
parser = argparse.ArgumentParser(description = "Get the ID abundance")
parser.add_argument('table', help="From blast table, outfmt 6 with this order: \
qaccver saccver ...etc - I get the first and the second column in order to obtain \
each gene with the ID or name of another table")
parser.add_argument('ID_table', help="Table which contains in the first column the \
same ID that the second column of the blast table and in the second column another ID")
parser.add_argument('cov', help="Coverage file with header")
parser.add_argument('out', help="Name/ID - coverage")
args = parser.parse_args()
# --------------------*********---------------------
# Read the Table
contigs = {}
with open(args.table, 'r') as file:
    for line in file:
        key = line.split("\t")[0]
        value = line.split("\t")[1]
        contigs[key] = value
IDs = {}
with open(args.ID_table, 'r') as file:
    for line in file:
        key = line.split("\t")[0]
        value = line.split("\n")[0].split("\t")[1]
        IDs[key] = value
# --------------------*********---------------------
# First I create a dictionary which contanis:
# ID or name AND BIN space separated
table = {}
for b in contigs:
    name = IDs[contigs[b]]
    if name not in table:
        table[name] = []
    bin = b.split("_k99")[0]
    table[name].append(bin)

# --------------------*********---------------------
# Read coverage table
coverage = pd.read_csv(args.cov, sep="\t", index_col=0)
# --------------------*********---------------------
# Get for every key in table2 the total abundance
dictionary = {}
for key in table:
    lista = []
    list_of_values = coverage.loc["Bin.100"] * 0
    for bin in table[key]:
        list_of_values = list_of_values +  coverage.loc[bin]
    for i in list_of_values:
        lista.append(i)
    dictionary[key] = lista
# --------------------*********---------------------
# Print the output
with open(args.out, 'w') as handle:
    #print('\n'.join(str(i) for i in dictionary), file=handle)
    handle.write(str(dictionary))
