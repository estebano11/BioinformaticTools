#!/usr/bin/env python3

### --- Import libraries ---- #
import argparse
import os
import numpy as np
from itertools import groupby
import pandas as pd
import cobra as cb
import micom as mc
from micom.community import Community
from joblib import Parallel, delayed
config = cb.Configuration()
config.solver = 'cplex'
### --- Import libraries ---- #

# Description of the script
parser = argparse.ArgumentParser(description = "Run micom in selected samples.\
    Warnings: - All the mags needs to have the same name for each table (abundance,fasta,models, etc)")
parser.add_argument('working_directory', type=str, default=".", 
                    help="Root path where all your directories and files are. Default: current directory")
parser.add_argument('path_fasta', help="Path where are the fasta mags/proteins used for generate the models")
parser.add_argument('path_models', help="path to models")
parser.add_argument('abundance', help="Abundance table with each MAG as row and samples in columns")
parser.add_argument('swap', type=int, default=0, 
                    help="If the samples in the abundance table are the rows and mags in columns set swap to 1. Default = 0")
parser.add_argument('simulations_out', help="Whole path to save the simulations")
parser.add_argument('fasta_extension', help="Extension used for fasta files (fa, fna, fasta, faa, etc)")
parser.add_argument('coverage_threshold', type=float, default=0.1, 
                    help="Threshold for the abundance, default=0.1")
parser.add_argument('check_growth', type=str, default="True", 
                    help="Check if the models grow without constraints, default=True")
parser.add_argument('num_threads', type=int, default=20, 
                    help="Number of cores to use for the simulations. Default=20")
args = parser.parse_args()

# --------------------*********--------------------- #
# --------------------*********--------------------- #
# --------------------*********--------------------- #
# --------------------*********--------------------- #
directory = os.path.abspath(args.working_directory)
dir_fasta = directory+"/"+args.path_fasta
dir_models = directory+"/"+args.path_models
swap = args.swap
if swap == 0:
    abundance = pd.read_csv(directory+"/"+args.abundance, sep="\t", index_col=0)
elif swap == 1:
    abundance = pd.read_csv(directory+"/"+args.abundance, sep="\t", index_col=0)
    abundance = abundance.T
else:
    print ("Swap is not well defined")
    print (swap)
    exit()
# check the extension
fasta_extension = args.fasta_extension
if "." not in fasta_extension:
    fasta_extension = "."+fasta_extension
coverage_threshold = args.coverage_threshold
check_growth = args.check_growth
simulations = directory+"/"+args.simulations_out

num_threads = args.num_threads
# --------------------*********--------------------- #
# --------------------*********--------------------- #
# --------------------*********--------------------- #
# --------------------*********--------------------- #
## ---- ##
# Define functions
def CreateMicomTable(abundance,dir_models,dir_fasta):
    """Create 'database_community' table for micom. The structure is the following:\
        id              fasta                     xml                   sample   abundance\
    0   mag_1   path/to/fasta/mag_1.fasta   path/to/models/mag_1.xml    sample 1    2.0"""
    
    database_community = pd.DataFrame(columns=["id","fasta","xml","sample","abundance"])
    for mag in abundance.index:
        fasta = dir_fasta+mag+fasta_extension
        xml = dir_models+mag+".xml"
        for sample in abundance.columns:
            cov = abundance.loc[mag][sample]
            database_community.loc[len(database_community)] = (mag,fasta,xml,sample,cov)
    return database_community

def CheckTableIntegrity(database_community,dir_models):
    """Check if the table provided has all the models"""
    missing = False
    missing_files = []
    models = list(database_community["xml"].unique())
    for xml in models:
        if os.path.isfile(xml):
            pass
        else:
            missing_files.append(xml)
            missing = True
    return missing,missing_files

def CreateTaxonomiesTable(database_community):
    """Create a list of dataframes to separate in samples"""
    database_community2 = database_community.rename(columns={"xml":"file"})
    database_community2.drop(columns=["fasta"], inplace=True)
    taxonomies = []
    samples = list(database_community2["sample"].unique())
    for sample in samples:
        tax = database_community2[database_community2["sample"] == sample]
        tax = tax[tax["abundance"] > coverage_threshold]
        taxonomies.append(tax)
    return taxonomies

    
def CheckGrowthModels(database_community):
    """Check if all the models grow without constraints"""
    models = list(database_community["xml"].unique())
    models_growth = {}
    all_growing = True
    for model in models:
        m = cb.io.read_sbml_model(model)
        s = m.optimize()
        gr = s.objective_value
        models_growth[model] = gr
        if gr == 0.0:
            all_growing = False
    
    return all_growing,models_growth

def runSimulations(tax):
    sample = list(tax["sample"])[0]
    print (f"Running sample {sample} with {len(tax)} species")
    com = Community(tax,solver="cplex")
    #com = Community(tax)
    com.to_pickle(f"{simulations}{sample}_{coverage_threshold}.pickle")
    return f"sample {sample} finished"
# --------------------*********--------------------- #
# --------------------*********--------------------- #
# --------------------*********--------------------- #
# --------------------*********--------------------- #
# Run script
# Create the table for micom and storage in the output of simulations
database_community = CreateMicomTable(abundance,dir_models,dir_fasta)
database_community.to_csv(f"{simulations}database_community.tsv",sep="\t")

# Check the integrity of the table
missing, missing_files = CheckTableIntegrity(database_community,dir_models)
if missing == True:
    # means that I have some files missing in the directory
    print ("Script aborted with errors, there are missing files")
    for file_ in missing_files:
        print (file_)
    exit()
else:
    print ("All files are present")
    
# Check if the models grow
if check_growth == "True":
    all_growing,models_growth = CheckGrowthModels(database_community)
    if all_growing == False:
        print ("Script aborted with errors")
        print ("Not all models are growing")
        for model in models_growth.keys():
            gr = models_growth[model]
            if gr == 0:
                print (model)
        exit()
    else:
        print ("All models are growing")
else:
    print ("I am not checking if all models are growing")

# create taxonomies dataframe
taxonomies = CreateTaxonomiesTable(database_community)
# run simulations
print ("Running simulations")
for tax in taxonomies:
	print (" --- ")
	runSimulations(tax)
#Parallel(n_jobs=num_threads)(delayed(runSimulations)(tax) for tax in taxonomies)
#runSimulations(taxonomies[0])
print ("Script finished")
