#!/usr/bin/env python3

import sys
import re
import csv

f = open(sys.argv[1])
line = f.readline()

while line != '':
	id=line
	line = f.readline()[:-1]
	if len(line) >= int(sys.argv[2]):
		print(id[:-1])
		print(line)
	line = f.readline()

f.close()
