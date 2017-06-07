#!/usr/bin/python

import sys

f = open(sys.argv[1])
line = f.readline()

while line[0] == "@":
	line = f.readline()

while line != '':
	t = line.split("\t");
	alLen = len(t[9])
	if alLen >= int(sys.argv[2]):
		print(line[:-1])
		#open(sys.argv[3] + 
	line = f.readline()

f.close()
