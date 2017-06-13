#!/usr/bin/python

import sys
import subprocess

f = open(sys.argv[1])

line = f.readline()
while line[0] == "@":
	line = f.readline()

finalString = ""
line = f.readline()
if line != "":
	t = line.split("\t");
	curFile = t[2]
	out = open(sys.argv[3] + curFile, "w")
	alLen = len(t[9])
	if alLen >= int(sys.argv[2]):
		finalString = finalString + line
	line = f.readline()	

while line != '':
	t = line.split("\t");
	alLen = len(t[9])
	if alLen >= int(sys.argv[2]):
		if curFile != t[2]:
			out.write(finalString)
			finalString = ""
			out.close()
			curFile = t[2]
			out = open(sys.argv[3] + curFile, "w")
		finalString = finalString + line
	line = f.readline()

out.write(finalString)
out.close()
f.close()
