#!/usr/bin/env python3
from Bio import Entrez
import argparse
import csv
# --------------------*********---------------------
# Description of the script
parser = argparse.ArgumentParser(description = "Get taxonomy information from id")
parser.add_argument('email', help="Entrez e-mail")
parser.add_argument('table', help="1 column table containing the id")
parser.add_argument('out', help="two column table: Id-Taxonomy")
args = parser.parse_args()
# --------------------*********---------------------
Entrez.email = args.email
# Read the table
IDs = {}
with open(args.table, 'r')as file:
    for f in file:
        key = str(f).split("\n")[0]
        print (key)
        handle = Entrez.efetch(db="nucleotide", id=key, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        value = record[0]['GBSeq_taxonomy']
        IDs[key] = value
w = csv.writer(open(args.out, 'w'))
for key, val in IDs.items():
        w.writerow([key,val])
