#!/usr/bin/env python
import os
import sys
import re
import argparse
import random
import collections




parser = argparse.ArgumentParser(description="""           
Description
-----------
    This script simulates single-end reads from the population genome""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Prerequisites
-------------
    python version 3+

Authors
-------
    Robert Kofler 
""")



parser.add_argument("--rm", type=str, required=True, dest="rm", default=None, help="the repeatmasker file")
args = parser.parse_args()

teh=collections.defaultdict(lambda:collections.defaultdict(lambda:0))
temax=collections.defaultdict(lambda:0)
speciesset=set([])



for line in open(args.rm):
	line=line.rstrip().lstrip()
	a=re.split("\s+",line)
	#print(a)
	#  0	1	2		3				4		5		6		7		8	9	10			11		12		13		14		
	# rm	294	16.87	contig_1005_1	47209	47291	C	Shellder	599	681	0.01	 A.communis
	# rm	234	15.01	contig_1005_1	47231	47332	C	Shellder	599	692	0.01	 A.communis
	score,te,species=int(a[1]),a[7],a[11]
	speciesset.add(species)
	if score> temax[te]:
		temax[te]=score
	if score > teh[te][species]:
		teh[te][species] =score

for te in temax.keys():
	dmelscore=temax[te]
	if dmelscore<1:
		continue
	for species in speciesset:
		relscore=float(teh[te][species])/float(dmelscore)
		print("{0}\t{1}\t{2}".format(te,species,relscore))
