'''
Author: RICCARDO PIANEZZA

This script is built to solve a common problem encountered with DeviaTE: TE sequences can be populated by
satellites, poly-A or other low complexity regions. This leads to problems with the mapping and usually to
an overestimation of the copynumber due to many reads mapping on these regions, as well as strechted plots.

To solve this, this script takes as input a folder containing the DeviaTE output files for a transposon for
all the samples, an interval of bases (the low complexity region, es. 23-145) and an output folder.

The script "flatten" the region to the average coverage of all the other positions in the TE in the sample,
re-calculates the estimated copynumber and HQ copynumber and write a new file with the new values.
It distributes the reads in the region equally among the 4 bases, so DO NOT USE THIS SCRIPT FOR SNP ANALYSIS
or simply remove the low complexity region for such goal.
'''

import os
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('input', metavar='input', type=str, help='path to input folder with all deviaTE outputs')
parser.add_argument('start', metavar='input', type=int, help='start of the low complexity region')
parser.add_argument('end', metavar='input', type=int, help='start of the low complexity region')
parser.add_argument('output', metavar='output', type=str, help='path to output folder')
args = parser.parse_args()

for filename in os.listdir(args.input):
    with open(args.input+str(filename), 'r') as input, open(args.output+str(filename), 'w') as output:
        cn_after = []
        cn_after_hq = []
        for line in input:
            if line.startswith("#"):
                continue
            else:
                base = int(line.split(' ')[2])
                if (base > args.start) and (base < args.end):
                    continue
                else:
                    coverage = line.split(' ')[8]
                    coverage_hq = line.split(' ')[10]
                    cn_after.append(float(coverage))
                    cn_after_hq.append(float(coverage_hq))
                    new_cn = sum(cn_after)/len(cn_after)
                    new_cn_hq = sum(cn_after_hq)/len(cn_after_hq)
        input.seek(0)
        for i, line in enumerate(input):
            if i == 0 or i == 2 or not ((i > args.start + 2 and i < args.end + 2) or i == 1):
                output.write(line)
            elif i == 1:
                new_line = "# insertions/haploid: " + str(new_cn) + " or " + str(new_cn_hq) + " (hq coverage only)"
                output.write(new_line + "\n")
            else:
                new_line = " ".join(line.split(' ')[:4]) + f" {new_cn/4} {new_cn/4} {new_cn/4} {new_cn/4} {new_cn}" + " " + line.split(' ')[9] + f" {new_cn_hq} {' '.join(line.split(' ')[11:])}"
                output.write(new_line)