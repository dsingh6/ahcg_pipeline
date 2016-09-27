#!/usr/bin/env python
import os
import sys

#usage: python master.py <bed> <vcf> <out.ext>

def bedvariants(bed,variants):	
	#parse bed file
	f=open(bed,"r")
	lines=f.readlines()
	f.close()
	d=[]
	for line in lines:
		if len(line)<3:
			continue
		splitline=line.split("\t")
		d.append([splitline[0],int(splitline[1]),int(splitline[2].strip())])
	#parse vcf file
	f=open(variants,"r")
	lines=f.readlines()
	f.close()
	varlist=[]
	for line in lines:
		if line.startswith("#"):
			continue
		splitline=line.split("\t")
		for item in d:
			if item[0] == splitline[0]:
				if int(splitline[1])> item[1] and int(splitline[1])<item[2]:
					varlist.append(line)
	varset=set(varlist)
	return varset

varset=bedvariants(sys.argv[1],sys.argv[2])
f=open(sys.argv[3],"w")
for var in varset:
	f.write(var)
f.close()
