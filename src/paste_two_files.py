#!/usr/bin/env python3
import argparse
from itertools import groupby
import pandas as pd

# --------------------*********---------------------
# Description of the script
parser = argparse.ArgumentParser(description = "Script to paste two tables (tab separated)")
parser.add_argument('table1', help="table1")
parser.add_argument('index1', help="name of the column index of table1")
parser.add_argument('table2', help="table2")
parser.add_argument('index2', help="name of the column index of table2")
parser.add_argument('out', help="coverage output")
args = parser.parse_args()

# --------------------*********---------------------
table1 = pd.read_csv(args.table1, sep="\t", index_col=args.index1)
del table1.index.name
table2 = pd.read_csv(args.table2, sep="\t", index_col=args.index2)
del table2.index.name

table_out = pd.concat([table1,table2], axis=1, join='inner')
table_out.to_csv(args.out, sep="\t",index_label="Bin")
