#!/usr/bin/python3

import sys
import re

f = open(sys.argv[1])

id = f.readline()[:-1]
while id != "":
	seq = f.readline()[:-1]
	splits = re.split('[acgtn]+', seq)
	i = 1
	for s in splits:
		print(id + "_" + str(i) + "\n" + s)
		i = i + 1
	id = f.readline()[:-1]
f.close()
