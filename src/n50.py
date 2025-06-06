#!/usr/bin/env python

# calculate N50 from fasta file
# N50 = contig length so that half of the contigs are longer and 1/2 of contigs are shorter

# Check the dependencies
dependencies = ["sys","os","itertools","numpy","argparse"]
# =============================================================== #
# import libraries
import sys
import os
from itertools import groupby
import numpy
import argparse
# =============================================================== #
parser = argparse.ArgumentParser(description="My Project - Run Specific Programs")
parser.add_argument("program", choices=["coverage", "aai", "depth"], help="Specify the program to run")
parser.add_argument("-i", "--input", help="Input filename")
parser.add_argument("-o", "--output", help="Output directory")
args = parser.parse_args()


# =============================================================== #
out_dir = args.out_dir
# check if the out_dir finish with slash
if out_dir[-1] != "/":
	out_dir = out_dir+"/"

# =============================================================== #

lengths = []

with open(sys.argv[1]) as fasta:
	# parse each sequence by header: groupby(data, key)
	faiter = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))
	GC = 0
	nombres = []
	GC_total = []
	for record in faiter:
        # join sequence lines to one.
		for line in record:
			if line[0] == ">":
				header1 = line.replace(">","")
				header = header1.replace('\n','')
		seq = "".join(s.strip() for s in faiter.next())
		#print header, seq
		gc = 0
		for nucleotide in seq:
			if nucleotide == "G" or nucleotide == "C" or nucleotide == "g" or nucleotide == "c":
				gc += 1
		gc_2 = (gc/float(len(seq)))*100
		GC = format(gc_2, '.2f')
		GC_total.append(float(GC))
		nombres.append(header + '\t' + str(len(seq)) + '\t' + str(GC))
		lengths.append(len(seq))

# N50
# the length of the shortest contig such that the sum of contigs of equal
# length or longer is at least 50% of the total length of all contigs

# sort contigs longest>shortest
all_len=sorted(lengths, reverse=True)
csum=numpy.cumsum(all_len)


Contigs = len(lengths)
GC_promedio = 0
for i in GC_total:
	GC_promedio = GC_promedio + i
gc_average = format(GC_promedio / Contigs, '.2f')
n2=int(sum(lengths)/2)

# get index for cumsum >= N/2
csumn2=min(csum[csum >= n2])
ind=numpy.where(csum == csumn2)
n50 = all_len[int(ind[0])]

Genome = str(sys.argv[1].split('/')[-1])
print (Genome)

# write lengths to file
with open(out_dir+Genome + '.stadistic', 'w') as handle:
	handle.write("Genome: " + Genome + '\t'"N50: %s" % n50 + '\t' "Number of contigs: %s" %Contigs + '\t' "Genome size: %s" %int(sum(lengths)) + '\t' "GC: %s" %gc_average)
	handle.write('\n'"Contig" + '\t' + "Length" + '\t' + "GC"'\n')
	handle.write('\n'.join(str(i) for i in nombres))
