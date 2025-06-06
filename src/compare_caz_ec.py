#!/usr/bin/env python3

import argparse
import os
import numpy as np
from itertools import groupby
import pandas as pd

# --------------------*********---------------------
# Description of the script
parser = argparse.ArgumentParser(description = "Script to get a table of cazymes \
for each mag and a final table.")
parser.add_argument('ec_gh', help="File with EC-GHs ")
parser.add_argument('caz', help="Output of the HMM of the MAG with cazy database")
parser.add_argument('eggnong', help="Output of the eggnong analysis")
parser.add_argument('out', help="prefix for output")
args = parser.parse_args()

# --------------------*********---------------------
# read table with EC - GHs
ec_table_tmp = pd.read_csv(args.ec_gh, sep="\t", names=["EC","GH"])
ec_table = pd.DataFrame(columns=["EC","GH"])

for i in ec_table_tmp.itertuples():
    ec = i[1]
    ghs = i[2]
    if " " in ghs:
        for gh in ghs.split(" "):
            ec_table = ec_table.append({"EC":ec,"GH":gh}, ignore_index=True)
    else:
        ec_table = ec_table.append({"EC":ec,"GH":ghs}, ignore_index=True)

# Read file of the results of dbCAN for the mag and the eggnog results
# cols = ["X","target name","accession","tlen","query name","accession2","qlen","E-value","score","bias","# of c-Evalue","i-Evalue","score2","bias2","from","to",
#         "from2","to2","from3","to3","acc description of target"]

caz = pd.read_csv(args.caz,sep="\t", header=None)
eggnog = pd.read_csv(args.eggnong, sep="\t")
mag_table = pd.DataFrame(columns=["CAZ","EC_EGGNOG","EC_CAZ_MATCH","COUNT"])

cazy_match = {}
cazy_not_match = {}

def caz_match(ec_eggnog, ecs, cazy_match,cazy_not_match):
    EC = ec_eggnog
    if ec_eggnog in ecs:
        EC_CAZ_MATCH = "YES"
        if cazy not in cazy_match.keys():
            cazy_match[cazy] = [cazy,EC,"YES",1]
        else:
            value = cazy_match[cazy][3]
            value += 1
            cazy_match[cazy] = [cazy,EC,"YES",value]
    else:
        EC_CAZ_MATCH = "NO"
        if cazy not in cazy_not_match.keys():
            cazy_not_match[cazy] = [cazy,EC,"NO",1]
        else:
            value = cazy_not_match[cazy][3]
            value += 1
            cazy_not_match[cazy] = [cazy,EC,"NO",value]
    return cazy_match,cazy_not_match


for gene in caz[0].unique(): # iterate over the unique genes
    caz_gene = caz[caz[0] == gene]
    # Choose the one with lower e-value:
    idx = caz_gene[6].sort_values().head(n=1).index

    cazy = caz_gene.loc[idx[0]][3].replace(".hmm","")
    if "_" in cazy:
        cazy = cazy.split("_")[0]
    # get all the values of that cazyme
    EC = "-"

    if "GH" in cazy:
        GH = cazy.replace("GH","")       
        ecs = ec_table[ec_table["GH"] == GH]["EC"].values
        ec_eggnog = eggnog[eggnog["#query"] == gene]["EC"].values
        ec_eggnog = ec_eggnog[0]
        if "," in ec_eggnog:
            for ec_egg in ec_eggnog.split(","):
                cazy_match,cazy_not_match = caz_match(ec_egg,ecs,cazy_match,cazy_not_match)
        else:
            cazy_match,cazy_not_match = caz_match(ec_eggnog,ecs,cazy_match,cazy_not_match)
        # if ec_eggnog in ecs:
        #     EC = ec_eggnog
        #     EC_CAZ_MATCH = "YES"
        #     if cazy not in cazy_match.keys():
        #         cazy_match[cazy] = [cazy,EC,"YES",1]
        #     else:
        #         value = cazy_match[cazy][3]
        #         value += 1
        #         cazy_match[cazy] = [cazy,EC,"YES",value]
        # else:
        #     EC = ec_eggnog
        #     EC_CAZ_MATCH = "NO"
        #     if cazy not in cazy_not_match.keys():
        #         cazy_not_match[cazy] = [cazy,EC,"NO",1]
        #     else:
        #         value = cazy_not_match[cazy][3]
        #         value += 1
        #         cazy_not_match[cazy] = [cazy,EC,"NO",value]

for key in cazy_match.keys():
    cazy, EC, ec_caz_match, count = cazy_match[key]
    mag_table = mag_table.append({"CAZ":cazy,"EC_EGGNOG":EC,"EC_CAZ_MATCH":ec_caz_match,"COUNT":count}, ignore_index=True)
for key in cazy_not_match.keys():
    cazy, EC, ec_caz_match, count = cazy_not_match[key]
    mag_table = mag_table.append({"CAZ":cazy,"EC_EGGNOG":EC,"EC_CAZ_MATCH":"NO","COUNT":count}, ignore_index=True)

mag_table.to_csv(args.out,sep="\t")
