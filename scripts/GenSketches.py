#!/usr/bin/env python3

import argparse as args
from Bio import SeqIO
from math import floor

#This script reads sequences from FASTA files and outputs their FracMinHash sketches

#Mapping between a nucleotide character and its bitwise representation
NT_IN_BITS = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

#This function calculates the hash of a bitwise k-mer representation. The function is influenced by the code of "The minimizer Jac-
#card estimator is biased and inconsistent." from Belbasi et al. (function "minimap2_hash(seed,v,mask)" in file 
#"minimap2_hash_uncompiled.py").
def getHash(kmer, mask):
	u = kmer & mask
	u = ((~u) + (u << 21)) & mask # u = (u<<21)-(u+1) = 77594587*u-1
	u = u ^ (u >> 24)
	u = ((u + (u << 3)) + (u << 8)) & mask # u *= 265
	u = u ^ (u >> 14)
	u = ((u + (u << 2)) + (u << 4)) & mask # u *= 21
	u = u ^ (u >> 28)
	u = (u + (u << 31)) & mask # u *= 2147483649

	return u

#This function calculates the FracMinHash sketch of a sequence. It is influenced by the code of "The minimizer Jaccard estimator is"
#biased and inconsistent.from Belbasi et al. (function "hash_sequence(seq,kmerSize,hashFunc,canonical=False)" in file 
#"jaccard_correction_test.py").
def calcSketch(seq, k, thres):
	sketch = []

	#Calculate mask
	mask = (4 ** k) - 1

	#Iterate of all k-mers in sequence
	for i in range(len(seq) - k + 1):
		kmerBits = 0

		#Get bit representation of k-mer
		for c in seq[i:i+k]:
			kmerBits = (kmerBits << 2) + NT_IN_BITS[c]

		#Calculate hash
		kmerHash = getHash(kmerBits, mask)

		#Add hash to sketch if it is small enough
		if kmerHash <= thres:
			sketch.append(kmerHash)

	return sketch

if __name__ == '__main__':
	#Setting up the argument parser
	parser = args.ArgumentParser(description="This script reads sequences from FASTA files and outputs their sketches.")
	parser.add_argument('-s', metavar='Seqs', type=str, nargs='+', required=True, help="Sequences in FASTA format")
	parser.add_argument('-k', metavar='KmerLen', type=int, default=15, help="Length of k-mers in sketch")
	parser.add_argument('-r', metavar='SampPar', type=float, default=0.1, help="Sampling parameter")

	arguments = parser.parse_args()
	minHashThres = floor(((4 ** arguments.k) - 1) * arguments.r)

	for s in arguments.s:
		for r in SeqIO.parse(open(s, 'r'), "fasta"):
			print(f">{r.description}")
			print(' '.join([str(h) for h in calcSketch(r.seq, arguments.k, minHashThres)]))
