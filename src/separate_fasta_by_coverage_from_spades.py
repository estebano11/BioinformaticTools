#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import sys
from itertools import groupby
import numpy as np
import re

with open(sys.argv[1]) as fasta:
with open("contigs.fasta") as fasta:
    faither = (x[1] for x in groupby(fasta, lambda line: line[0] == ">"))
    cov = []
    nombres = []
    for record in faither:
        for line in record:
            if line[0]==">":
                header1 = line.replace(">","")
                header = header1.replace('\n','')
                cov1 = re.sub(r".*cov_", r"", header)
        seq = "".join(s.strip() for s in faither.next())
        nombres.append(header + '\t' + str(len(seq)) + '\t' + str(cov1) + '\t' + str(seq))
        cov.append(cov1)


a = np.asarray(cov, dtype=np.float) #Tengo que pasar la lista a un array y convertir
                                   # la lista a float, ya que esta como strings
media = np.mean(a)                  # Calculo la media

# Aquellas secuencias que tengan un valor de cobertura +- 5 de la media
# me los quedo
secuencias = []
for item in nombres:
    if media - 5 < float(item.strip().split('\t')[2]) < media+5:
        secuencias.append (">" + item.strip().split('\t')[0] + '\n' + item.strip().split('\t')[3])

Bin = str(sys.argv[1].split('/')[-1])
with open(Bin + '.refinado', 'w') as handle:
    handle.write('\n'.join(str(i) for i in secuencias))

with open('array','w') as handle:
    for i in a:1
        handle.write(str(a))



