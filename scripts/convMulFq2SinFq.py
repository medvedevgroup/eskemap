#!/usr/bin/env python3

import sys
from Bio import SeqIO

#This scripts converts a multiple FASTQ file into single FASTQ files

outFilePrefix = sys.argv[1].split(".fastq")[0]
counter = 0

for r in SeqIO.parse(open(sys.argv[1], 'r'), "fastq"):
	SeqIO.write(r, open(f"{outFilePrefix}_{counter}.fastq", 'w'), "fastq")
	counter += 1
