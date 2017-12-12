#!/usr/bin/python3

import sys

f = open(sys.argv[1])
line = f.readline()

while line != '':
    id = line[:-1]
    line = f.readline()
    length = len(line) - 1

    s = id+"_"+str(length)
    print(s)
    print(line[:-1])

    line = f.readline()

f.close()
