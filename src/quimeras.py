# Program to count read per contig and MAG and detect if any read is shared with other MAG
import pandas as pd
import numpy as np
import argparse
import time


# ---------------------------------------------------- #
# Script's description
parser = argparse.ArgumentParser(description="Get the number of reads per contig and MAG.\n \
    Build 4 matrix: \n \
        - MAG + number of reads + reads aligned.\n\
        - MAG - number of reads shared with other MAGs.\n\
        - MAG - names of reads aligned with other MAGs.\n\
        Caution: The contigs are formated so the name of each one are: MAG-Name_Contig-Name.\n\
            For example: MAGs.10_NODE_592_length_71555_cov_65.571371 --> MAGs.10 --> MAG_Name and NODE... --> contig name")
parser.add_argument("sam", help='Sam formated file: \n\
    grep -v "^@" file.sam | cut -f1,3 | grep -v "\*" > file.txt")')
parser.add_argument("prefix", help="Prefix for the output")
args = parser.parse_args()
# ---------------------------------------------------- #

print ("-----------------------------------")
print ("       START       ")
start_time = time.time()

print ("START READING FILES")
df = pd.read_csv(args.sam, sep="\t", names=["read","mag"])
print ("FINISHED READING FILES")

print ("START FORMATING FILES")
df["mag"] = df["mag"].str.split("_NODE").str[0] # Replace the names of the contigs just to keep the mag name
df["mag"] = df["mag"].str.split("_k99").str[0]

df2 = df.groupby('mag').read.agg([('count', 'count'), ('read', ','.join)]) # Join the dataframe by mag names and count
print ("FINISHED FOMATING FILES")

mags = list(df2.index)
# df_number_of_reads = pd.DataFrame(index=args.MAGs)
df_reads_aligned = pd.DataFrame(index=mags)
df_number_of_reads_shared = pd.DataFrame(index=mags, columns=mags)
df_shared_reads = pd.DataFrame(index=mags, columns=mags)

d_reads = dict((el,0) for el in mags)
d_reads_aligned = dict((el,[]) for el in mags)

print ("START COUNTING")

mags_list = list(df2.index) # save a list of mags, iterate over the dataframe and delete the mag in the list to not repeat
total = len(mags_list)
for i in df2.itertuples():
    print ("{} % COMPLETED".format((total - len(mags_list))*100/total))
    mag = i[0]
    reads = set(i[2].split(","))
    mags_list.remove(mag)
    for mag2 in mags_list:
        reads2 = set(df2.loc[mag2]["read"].split(","))
        union = reads & reads2
        df_shared_reads.loc[mag,mag2] = list(union)
        df_shared_reads.loc[mag2,mag] = list(union)
        df_number_of_reads_shared.loc[mag2,mag] = len(union)
        df_number_of_reads_shared.loc[mag,mag2] = len(union)



# ---------------------------------------------------- #
# OUTPUT
prefix = str(args.prefix)
df_shared_reads.to_csv(prefix+"_shared_reads.tab", sep="\t")
df_number_of_reads_shared.to_csv(prefix+"_number_of_shared_reads.tab", sep="\t")
# df_number_of_reads.to_csv(prefix+"_number_of_reads", sep="\t")
df2.to_csv(prefix+"_reads_aligned.tab", sep="\t")

# ---------------------------------------------------- #
print ("-----------------/------------")
print("FINISHED IN %s seconds" % (time.time() - start_time))