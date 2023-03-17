#!/usr/bin/env python3

import argparse as args
from CalcLocInterSim import calcSketch, getHash, NT_IN_BITS
from MutateSeq import mutateSeq, NUCL_ALPHABET
from random import randrange, choice, seed
from sys import maxsize, stderr
from math import floor
from numpy import unique
from collections import deque
from Bio.Seq import Seq

#This script randomly generates a pair of sequences where one has been mutated according to a mutation process from the other. 
#Afterwards a sequencing process was simulated by introducing some additional errors. Sketches are calculated for both sequences and
#outputted.

#This function calculates the minimizer sketch of a sequence. It is influenced by the code of "The minimizer Jaccard estimator is 
#biased and inconsistent." from Belbasi et al. (function "winnowed_minimizers_linear(perm,windowSize)" in file 
#"winnowed_minimizers.py").
def calcMiniSketch(seq, k, w):
	sketch = []
	#A deque to store k-mer hashes inside the current window
	windowKmers = deque()
	mask = (4 ** k) - 1
	lastIdx = -1

	for i in range(len(seq) - k + 1):
		kmerBits = 0
		kmerBitsRevComp = 0
		windowBorder = i - (w - 1)

		#Get bit representation of k-mer
		for c in seq[i:i+k]:
			kmerBits = (kmerBits << 2) + NT_IN_BITS[c]

		#Get bit representation of k-mer's reverse complement
		for c in str(Seq(seq[i:i+k]).reverse_complement()):
			kmerBitsRevComp = (kmerBitsRevComp << 2) + NT_IN_BITS[c]

		#If a k-mer is its own reverse complement we skip it
		if kmerBits == kmerBitsRevComp:
			continue

		#Depending on which hash is smaller we consider either a k-mer or its reverse complement per position
		if kmerBits < kmerBitsRevComp:
			#A k-mer is a pair of k-mer's start position and its hash
			kmer = (i, getHash(kmerBits, mask))
		else:
			#A k-mer is a pair of k-mer's start position and its hash
			kmer = (i, getHash(kmerBitsRevComp, mask))

		#Remove all k-mers with a hash value larger than the newly calculated one
		while (len(windowKmers) > 0) and (windowKmers[-1][1] > kmer[1]):
			windowKmers.pop()

		#Save new k-mer as window k-mer
		windowKmers.append(kmer)

		#Remove k-mer if it is not any longer inside the window
		while (len(windowKmers) > 0) and (windowKmers[0][0] < windowBorder):
			windowKmers.popleft()

		#As soon as we have seen a first full window of k-mers choose a minimizer
		if (windowBorder >= 0) and (len(windowKmers) > 0):
			#We do not choose the same minimizer for a second time
			if lastIdx != windowKmers[0][0]:
				sketch.append(windowKmers[0][1])
				lastIdx = windowKmers[0][0]

			while len(windowKmers) > 1 and windowKmers[0][1] == windowKmers[1][1]:
				windowKmers.popleft()
				sketch.append(windowKmers[0][1])
				lastIdx = windowKmers[0][0]

	#If our sequence was too small to get a full window of k-mers to consider take the smallest one found so far
	if windowBorder < 0 and len(windowKmers) > 0:
		sketch.append(windowKmers[0][1])

		while len(windowKmers) > 1 and windowKmers[0][1] == windowKmers[1][1]:
			windowKmers.popleft()
			sketch.append(windowKmers[0][1])

	return sketch

if __name__ == '__main__':
	#Testing
	from Bio import SeqIO
	refSeq = str([r for r in SeqIO.parse(open("../../simulations/genomes/t2thumanChrY.fasta", 'r'), "fasta")][0].seq)
	blks = {}#{int(l): None for l in open("../highAbundKmersMiniK15w10Lrgr100BtStrnds.txt", 'r')}
	for h in [k for k in calcMiniSketch(refSeq, 15, 10) if not k in blks]:
		print(h)
	exit(0)

	#Setting up the argument parser
	parser = args.ArgumentParser(description="This script generates sketches of a sequence and its mutated copy in .sk format.")
	parser.add_argument('-l', metavar='SeqLen', type=int, required=True, help="Length of sequences")
	parser.add_argument('-m', metavar='SubRate', type=float, required=True, help="Substitution rate")
	parser.add_argument('-d', metavar='DelRate', type=float, required=True, help="Deletion rate")
	parser.add_argument('-i', metavar='InsLen', type=float, required=True, help="Mean insertion length")
	parser.add_argument('-se', metavar='SubErr', type=float, required=True, help="Substitution error")
	parser.add_argument('-de', metavar='DelErr', type=float, required=True, help="Deletion error")
	parser.add_argument('-ie', metavar='InsErr', type=float, required=True, help="Mean insertion length error")
	parser.add_argument('-s', metavar='Seed', type=int, default=randrange(maxsize), help="Random seed to use")
	parser.add_argument('-k', metavar='KmerLen', type=int, default=9, help="The k-mer length to use")
	parser.add_argument('-r', metavar='HashRatio', type=float, default=0.1, help="The ratio of hash values from the set of all " + \
		"possible hash values to be included into a sketch")
	parser.add_argument('-H', metavar='HashMode', type=str, default="FracMin", help="Hashing method to be used for sketches")
	parser.add_argument('-w', metavar='WindowSize', type=int, default=10, help="Window size for minimizer sketching approach")
	parser.add_argument('-b', metavar='Blacklist', type=args.FileType('r'), required=False, help="File containing blacklisted " + \
		"k-mers not to appear inside the sketch")

	arguments = parser.parse_args()

	#Check if rates and insertion length were chosen reasonably
	if arguments.m < 0 or arguments.m > 1:
		print("ERROR: Substitution rate should be between 0 and 1", file=stderr)
		exit(-1)

	if arguments.d < 0 or arguments.d > 1:
		print("ERROR: Deletion rate should be between 0 and 1", file=stderr)
		exit(-1)

	if arguments.i < 0:
		print("ERROR: Mean insertion length should be positive", file=stderr)
		exit(-1)

	#Generate sequence whose sketch has no duplicates
	seed(arguments.s)
	minHashThres = floor(((4 ** arguments.k) - 1) * arguments.r)

	seq = ""

	for i in range(arguments.l):
		seq += choice(NUCL_ALPHABET)

	if arguments.H == "FracMin":
		sketch = calcSketch(seq, arguments.k, minHashThres)
	elif arguments.H == "Mini":
		sketch = calcMiniSketch(seq, arguments.k, arguments.w)
	else:
		print("ERROR: Unrecognized sketching mode", file=stderr)
		exit(-1)

	#Check if blacklist was provided and apply filtering
	if arguments.b:
		blKmers = {int(l): None for l in arguments.b}
		sketch = [k for k in sketch if not k in blKmers]

		#Testing
		print(sketch)

	#Mutate sequence
	mutSeq = mutateSeq(seq, arguments.m, arguments.d, arguments.i)
	seqSeq = mutateSeq(mutSeq, arguments.se, arguments.de, arguments.ie)

	if arguments.H == "FracMin":
		mutSeqSketch = calcSketch(seqSeq, arguments.k, minHashThres)
	else:
		mutSeqSketch = calcMiniSketch(seqSeq, arguments.k, arguments.w)

	#Output sketches
	print(f">SketchPair_{arguments.s} original template")
	print(' '.join([str(h) for h in sketch]))
	print(f">SketchPair_{arguments.s} mutated and sequenced template")
	print(' '.join([str(h) for h in mutSeqSketch]))
