#!/usr/bin/python3

import sys

f = open(sys.argv[1])

id = f.readline()
while id != "":
	seq = f.readline()[:-1]
	print(id + seq.strip('acgtn'))
	id = f.readline()
f.close()
