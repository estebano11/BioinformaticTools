#!/usr/bin/env python3
import requests
import sys
import argparse
import pandas as pd
import csv
# --------------------*********---------------------
# Script's description
parser = argparse.ArgumentParser(description = "From blast table, class-names, class-go:id table  \
and coverage table, get the merop-abundance by class")
parser.add_argument('table', help="From blast table, outfmt with this order: \
contig merop - I get the first and the third column in order to obtain \
each gene with the ID or name of another table")
parser.add_argument('family', help="Two columns table: mer-family")
parser.add_argument('cov', help="Coverage file with header")
parser.add_argument('out', help="caz-coverage")
args = parser.parse_args()
# --------------------*********---------------------
# Read the Table
conversion = {'Aspartic':'A',
              'Cysteine':'C',
              'Glutamic':'G',
              'Metallo':'M',
              'Asparagine':'N',
              'Mixed': 'P',
              'Serine':'S',
              'Threonine':'T',
              'Unknown':'U',
              'Compund Peptidase':'X',
              'Inhibitors':'I'}
contigs = {}
with open(args.table, 'r') as file:
    for line in file:
        #key = line.split("\t")[0]
        key = (line.split("\t")[1])
        #value = (line.split("\t")[2]).split("_k99")[0]
        value = (line.split("\n")[0]).split("\t")[0].split("_k99")[0]
        if key not in contigs:
            contigs[key] = [value]
        else:
            contigs[key].append(value)
family = {}
with open(args.family, 'r') as file:
    for line in file:
        key = (line.split("\t")[0])
        value = (line.split("\n")[0]).split("\t")[1]
        family[key] = value

# GO = {}
# with open(args.go, 'r')as file:
#     for line in file:
#         key = line.split("\t")[1]
#         value = line.split("\t")[2]
#         GO[key] = value

# --------------------*********---------------------
# First I create a dictionary which contanis:
# ID or name AND BIN space separated
table = {}
for mer in contigs:
    s = family[mer][0]
    if s not in table:
        table[s] = []
        for bin in contigs[mer]:
            table[s].append(bin)
    else:
        for bin in contigs[mer]:
            table[s].append(bin)
new = {}
for key in table:
    for n in conversion:
        if conversion[n]==key:
            new_key = n
            new[n]=table[key]
table = new
# --------------------*********---------------------
# Read coverage table
coverage = pd.read_csv(args.cov, sep="\t", index_col=0)
# --------------------*********---------------------
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
