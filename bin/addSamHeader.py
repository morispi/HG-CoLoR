#!/usr/bin/env python3

import sys
import subprocess

f = open(sys.argv[1])

header = f.readline()[:-1]
while header != "":
	seq = f.readline()
	print("@SQ" + "\t" + header + "\t" + "LN:" + str(len(seq)))
	header = f.readline()[:-1]
print("@PG	ID:bwa	PN:bwa	VN:0.7.12-r1044	CL:bwa mem -t 16 /mnt/nfs/NGSData/referenceGenomes/Yeast/YEAST.fasta corYeast.split.fasta")
f.close()
