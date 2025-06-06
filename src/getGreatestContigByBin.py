#!/usr/bin/env python
import sys
from Bio import SeqIO
from itertools import groupby


archivo = str(sys.argv[1])
with open(sys.argv[1]) as fasta:
    faiter = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))
    largest_contig = ["none", 0]
    for record in faiter:
        # join sequence lines to one.
        for line in record:
            if line[0] == ">":
                header1 = line.replace(">","")
                header = header1.replace('\n','')
        seq = "".join(s.strip() for s in faiter.next())
        if len(seq) > largest_contig[1]:
            largest_contig = [header, len(seq)]
print archivo, "\t".join([str(i) for i in largest_contig])
