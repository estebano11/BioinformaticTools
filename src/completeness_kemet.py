#!/usr/bin/env python3
# --------------------*********---------------------
import argparse
import os
import pandas as pd
# Description of the script
parser = argparse.ArgumentParser(description = "Script for get Module completeness from KEMET output")
parser.add_argument('kemet', help="Full path of the kemet tsv files")
parser.add_argument('output', help="Name for the output file")
parser.add_argument('prefix', default="reportKMC_", help="Prefix of the kemet files, default: reportKMC_")
args = parser.parse_args()
# --------------------*********---------------------
out = args.output
kemet = args.kemet

module_completeness = pd.DataFrame()

for archive in os.listdir(kemet):
    # Check if is a tsv file:
    if "tsv" in archive:
        mag = archive.replace("reportKMC_","").replace(".tsv","")
        # create dataframe:
        columns = ["module","name","completeness","present_total","KOs present","KOs total"]
        df = pd.read_csv(kemet+archive, sep="\t",index_col=0, names=columns)
        for mod in df.index:
            present,total = df.loc[mod]["present_total"].split("__")
            if total > 3:
                completeness = total/present * 100
            elif present == total:
                completeness = 100
            elif present > 0:
                completeness = total/present * 100
            else:
                pass
            module_completeness.loc[mod][mag] = completeness
module_completeness.to_csv(out,sep="\t")
        