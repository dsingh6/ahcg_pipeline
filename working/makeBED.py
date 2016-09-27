#!/usr/bin/env python

import sys

#input refgene file
my_file = sys.argv[1]

#list of targets
target_list=['NM_032043',
'NM_007294',
'NM_000059',
'NM_004675',
'NM_001005862',
'NM_001080124',
'NM_000660',
'NM_000249',
'NM_000251',
'NM_000179',
'NM_000535',
'NM_002354',
'NM_000546',
'NM_000314',
'NM_000455',
'NM_004360',
'NM_024675',
'NM_001005735',
'NM_000044',
'NM_000051',
'NM_002485',
'NM_000465',
'NM_032043',
'NM_005732',
'NM_001164269',
'NM_058216',
'NM_002878']

#buffer to add/subtract to end/start of cds
start_buffer = 200
exon_buffer = 20
coord_list = []
init_start = 0
current_start = 0
current_end = 0
text = ""
revised_text = ""
i = 1
q = 1


firstline = True
with open(my_file) as f:
	the_text = f.readlines()

#iterate through target list
for target in target_list:
	for line in the_text:
		if target in line:
			props = line.split('\t')
			chrom = props[2]
			cds_start = props[6]
			cds_end = props[7]
			strand = props[3]
			coord_1 = props[9].strip()
			coord_2 = props[10].strip()
			all_coords = coord_1 + coord_2
				
	cds_s = int(cds_start)
	cds_e = int(cds_end)
	
	for coord in sorted(all_coords.split(',')):
		if coord.strip() == '':
			continue
		if i == 1:
			current_start = coord
			i = 2
			continue
		if i == 2:	
			current_end = coord
			i = 1
			text += chrom + "\t" + str(int(current_start)-exon_buffer) + "\t" + str(int(current_end)+exon_buffer) + "\n"
			continue
	
	#grab cds start and stop and replace start and end of gene coordinates with respective
	for line in text.split("\n"):
		if line.strip() != "":
			start = str(int(line.split("\t")[1])-start_buffer)
			stop = str(int(line.split("\t")[2])+start_buffer)
			chr = line.split("\t")[0]
			if cds_s > int(start) and cds_s < int(stop):
				revised_text += chr + "\t" + str(cds_s) + "\t" + stop + "\n"
				continue
			if cds_e > int(start) and cds_e < int(stop):
				revised_text += chr + "\t" + start + "\t" + str(cds_e) + "\n"
				break
			else:
				revised_text += chr + "\t" + start + "\t" + stop + "\n"
		
print revised_text
