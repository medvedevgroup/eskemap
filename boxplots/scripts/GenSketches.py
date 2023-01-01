#!/usr/bin/env python3

from sys import argv
from Bio import SeqIO
from math import floor
from CalcLocInterSim import calcSketch
from os.path import basename

#This script generates a sketch pair from a given template sequence and a read file

if __name__ == '__main__':
	K = 15
	R = 0.1
	MIN_HASH_THRES = floor(((4 ** K) - 1) * R)
	genome = basename(argv[1]).split('.fasta')[0]
	seed = basename(argv[2]).split('_s')[1].split('_d')[0]

	#Iterate over input files
	for i in range(1, 3):
		for r in SeqIO.parse(open(argv[i], 'r'), "fasta"):
			sketch = calcSketch(r.seq, K, MIN_HASH_THRES)

			if i == 1:
				print(f">ReadSketchPair_{genome}_{seed} original template")
			else:
				print(f">ReadSketchPair_{genome}_{seed} sequenced template")

			print(' '.join([str(h) for h in sketch]))
