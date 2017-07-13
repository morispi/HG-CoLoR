#!/usr/bin/python

import sys
import re
import csv

reads = open(sys.argv[1])
counts = open(sys.argv[2])
th = int(sys.argv[3])

id = reads.readline()[:-1]
while id != '':
	seq = reads.readline()[:-1]
	counts.readline()
	t = counts.readline().split(" ")
	i = 0
	bad = False
	while i < len(t) and not bad:
		if int(t[i]) < th:
			bad = True
		i = i + 1
	if not bad:
		print(id)
		print(seq)
	id = reads.readline()[:-1]

reads.close()
counts.close()

