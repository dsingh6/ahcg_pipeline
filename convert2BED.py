#!/usr/bin/env python

from re import *
import os
import sys

#filename = sys.argv[1]
#outfile = sys.argv[2]

f = open(filename, "r")
fo = open(outfile, "w")
a = 1
for line in f:
	splitLine = line.split()
	chrm = splitLine[2]
	strand = splitLine[3]
	name = splitLine[1]
	starts = splitLine[9]
	startSplit = starts.split(',') 
	ends = splitLine[10]
	endSplit = ends.split(',')
	stop = len(startSplit) - 1
	for i in range(0,stop):
		fo.write(chrm + "\t"  + startSplit[i] + "\t" + endSplit[i] + "\t" + name + "\t" + "." + "\t" + strand + "\n")	
