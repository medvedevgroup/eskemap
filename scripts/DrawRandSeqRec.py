#!/usr/bin/env python3

import sys
from Bio import SeqIO
from random import randint

#This script draws a random sequence from a fastq file and saves it in FASTA format

seqRecs = [r for r in SeqIO.parse(open(sys.argv[1], 'r'), "fastq")]
draw = randint(0, len(seqRecs) - 1)
SeqIO.write(seqRecs[draw], open(sys.argv[2], 'w'), "fasta")
