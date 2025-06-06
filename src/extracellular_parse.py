import os
import argparse
import pandas as pd
import re

######################################################
# Description of the script
parser = argparse.ArgumentParser(description = "Script to get the protein envelope from psortb")
parser.add_argument('file', help="psortb output file")
parser.add_argument('gene_mags', help="File with gene Id + mag")
parser.add_argument('out', help="output name")
args = parser.parse_args()
######################################################
file = args.file
gene_mags = args.gene_mags
gene_mags_df = pd.read_csv(gene_mags, index_col=0, names=["mag"], sep="\t")

df = pd.DataFrame(columns=["mag", "envelope", "gene"])

next = False
with open (file) as f:
    for line in f:
        # line = re.sub('\s+',' ',line)
        if next == True:
            envelope = re.sub('\s+',' ',line).replace("\n","")
            df = df.append({"mag":mag, "envelope":envelope, "gene":gene},ignore_index=True)
            next = False
            continue
        if line.startswith("SeqID"):
            gene = line.replace("\n","").replace("SeqID: ","")
            mag = gene_mags_df.loc[gene]["mag"]
        else:
            if "Final Prediction" in line:
                next = True

df.to_csv(args.out,index=False, sep="\t")
