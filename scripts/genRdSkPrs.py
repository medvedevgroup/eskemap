#!/usr/bin/env python3

from sys import argv, stderr
from Bio import SeqIO
from math import floor, ceil
from CalcLocInterSim import calcSketch
from os.path import basename
import gzip

#This script generates a sketch pair from a given template sequence and a read file

#This function parses a FASTA file and returns its content
def readFasta(path):
	return [r for r in SeqIO.parse(open(path, 'r'), "fasta")]

if __name__ == '__main__':
	K = 15
	R = 0.1
	MIN_HASH_THRES = floor(((4 ** K) - 1) * R)
	genome = basename(argv[2]).split('_c')[0]
	seed = basename(argv[2]).split('_s')[1].split('_d')[0]
	refSt = -1

	for l in gzip.open(argv[2], 'rt'):
		if l.startswith("s ref"):
			for e in l.split(' '):
				if e.isnumeric():
					if refSt >= 0:
						refL = int(e)
						break
					else:
						refSt = int(e)
		elif l.startswith("s S"):
			rdSeq = l.strip().split(' ')[-1].replace('-', '')
			break

	#Testing
	# print("refL:", refL)
	# print("refSt:", refSt)

	refSt = max(0, refSt - ceil((len(rdSeq) - refL) / 2))
	oRef = readFasta(argv[1])[0].seq[refSt:refSt + len(rdSeq)]

	#Sanity check: Outputted reference and read sequences short be of same length
	if len(oRef) != len(rdSeq):
		print("ERROR: Reference and read sequence do not have the same length!", file=stderr)
		exit(-1)

	#Testing
	# print(len(oRef))
	# print(len(rdSeq))
	# exit(0)

	print(f">ReadSketchPair_{genome}_{seed} original template")
	print(' '.join([str(h) for h in calcSketch(oRef, K, MIN_HASH_THRES)]))
	print(f">ReadSketchPair_{genome}_{seed} sequenced template")
	print(' '.join([str(h) for h in calcSketch(rdSeq, K, MIN_HASH_THRES)]))

# 	#Iterate over input files
# 	for i in range(1, 3):
# 		for r in SeqIO.parse(open(argv[i], 'r'), "fasta"):
# 			sketch = calcSketch(r.seq, K, MIN_HASH_THRES)

# 			if i == 1:
# 				print(f">ReadSketchPair_{genome}_{seed} original template")
# 			else:
# 				print(f">ReadSketchPair_{genome}_{seed} sequenced template")

# 			print(' '.join([str(h) for h in sketch]))

# import gzip
# from math import ceil

# refSt = -1
    
# for l in gzip.open("../simulations/reads/randSeq_l110_rid0_cP6C4_ep6:50:54_s7361077429744071834_d1_l100-100.maf.gz"\
#                    , 'rt'):
#     l = l.strip()
    
#     if l.startswith("s ref"):
#         for e in l.split(' '):
#             if e.isnumeric():
#                 if refSt >= 0:
#                     refL = int(e)
#                     break
#                 else:
#                     print(e)
#                     refSt = int(e)
#     elif l.startswith("s S"):
#         rdSeq = l.split(' ')[-1].replace('-', '')
#         break

# print("RefSeq:", refSeq)
# print("len:", len(refSeq))
# print("RdSeq:", rdSeq)
# print("len:", len(rdSeq))
# print("refSt:", refSt)
# print("refL:", refL)
# cutStart = refSt - ceil((len(rdSeq) - refL) / 2)

# origRef = readFasta("../simulations/genomes/randSeq_l110_rid0.fasta")[0].seq[cutStart:cutStart + len(rdSeq)]

# print(origRef)
# print(len(origRef))