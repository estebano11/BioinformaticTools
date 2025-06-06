#!/usr/bin/env python3
# --------------------*********---------------------
import argparse
import os
import pandas as pd
# Description of the script
parser = argparse.ArgumentParser(description = "Script to join gene and contig names in htseq output")
parser.add_argument('folder', help="Full path of the tsv files folder - Avoid another tsv files in the same folder")
parser.add_argument('output_extension',default=".tab", help="Extension for the output, default: .tab")
parser.add_argument('add_one',default=False, help="Add 1 to all the values - True or False, default: False")
args = parser.parse_args()
# --------------------*********---------------------
folder = args.folder
extension = args.output_extension
add = args.add_one

for file_ in os.listdir(folder):
    if file_.endswith(".tsv"):
        lines = []
        with open(folder+file_)as f:
            for line in f:
                gene,contig,value = line.split("\t")
                value = float(value.replace("\n",""))
                g = gene.split("_")[1]
                contig_gene = contig+"_"+g
                if add == "True":
                    value += 1
                lines.append(contig_gene+"\t"+str(value)+"\n")
                
        output = file_.replace("tsv",extension)
        with open(output,'w')as out:
            for l in lines:
                out.write(l)