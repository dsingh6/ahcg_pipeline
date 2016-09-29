#!/usr/bin/env python

import sys
import getopt

inputfile = ''
outputfile = ''
opts, args = getopt.getopt(sys.argv[1:],"i:o:")

for opt, arg in opts:
	if opt == '-i':
		inputfile = arg
	elif opt == '-o':
  	        outputfile = arg
	else:
   	        print "Usage: ./BEDmaker.py -i filename -o filename"


o = open(outputfile, 'w')

with open(inputfile) as f:
	for line in f:
		gene = line
		aGene = gene.split()

		start = aGene[9]
		end = aGene[10]

		aStart = start.split(',')
		aEnd = end.split(',')


		for l in range(len(aStart) - 1):
					
			o.write(aGene[2] + "\t" + str(int(aStart[l]) - 20) + "\t" + str(int(aEnd[l]) + 20) + "\t" + aGene[1] + ':' + aGene[12] + "\t" + "NA" + "\t" + aGene[3])
			o.write("\n") 


f.close()
o.close()
