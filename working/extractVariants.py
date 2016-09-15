#!/usr/bin/env python

from re import *
import os
import sys

inFile = sys.argv[1]
#outFile = sys.argv[2]

# Open files to be read and written
f=open(inFile, "r")
#fo=open(outFile, "w")

for line in f:
	if not line.startswith('#'):
		splitLine = line.split()
		chrom = splitLine[0]
		start = splitLine[1]
		if chrom == '17' or chrom == 'chr17':
			start = int(start)
			if start >= 41196311 and start <= 41197819:
				print line
			elif start >= 41199659 and start <= 41199720:
				print line
			elif start >= 41201137 and start <= 41201211:
				print line
			elif start >= 41203079 and start <= 41203134:
				print line
			elif start >= 41209068 and start <= 41209152:
				print line
			elif start >= 41215349 and start <= 41215390:
				print line
			elif start >= 41215890 and start <= 41215968:
				print line
			elif start >= 41219624 and start <= 41219712:
				print line
			elif start >= 41222944 and start <= 41223255:
				print line
			elif start >= 41226347 and start <= 41226538:
				print line
			elif start >= 41228504 and start <= 41228631:
				print line
			elif start >= 41234420 and start <= 41234592:
				print line
			elif start >= 41242960 and start <= 41243049:
				print line
			elif start >= 41243451 and start <= 41246877:
				print line
			elif start >= 41247862 and start <= 41247939:
				print line	
			elif start >= 41249260 and start <= 41249306:
				print line
			elif start >= 41251791 and start <= 41251897:
				print line	
			elif start >= 41256138 and start <= 41256278:
				print line	
			elif start >= 41256884 and start <= 41256973:
				print line
			elif start >= 41258472 and start <= 41258550:
				print line
			elif start >= 41267742 and start <= 41267796:
				print line
			elif start >= 41276033 and start <= 41276132:
				print line
			elif start >= 41277287 and start <= 41277500:
				print line
