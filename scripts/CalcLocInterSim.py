#!/usr/bin/env python3

import argparse as args
from math import floor

#This script is supposed to find locally similar blocks inside sequences by considering their sketches and the so called "local intersection 
#similarity" measure

#Mapping between a nucleotide character and its bitwise representation
NT_IN_BITS = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

#This function calculates the hash of a bitwise k-mer representation. The function is influenced by the code of "The minimizer Jaccard estimator is 
#biased and inconsistent." from Belbasi et al. (function "minimap2_hash(seed,v,mask)" in file "minimap2_hash_uncompiled.py").
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

#This function calculates the sketch of a sequence. It is influenced by the code of "The minimizer Jaccard estimator is biased and inconsistent."
#from Belbasi et al. (function "hash_sequence(seq,kmerSize,hashFunc,canonical=False)" in file "jaccard_correction_test.py").
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

#This function calculates the "intersection similarity" for two intervals of two given strings
def calcSim(skA, s, e, skB, t, f):
	setA = set(skA[s:e+1])
	setB = set(skB[t:f+1])

	return len(setA.intersection(setB)) - len(setA.difference(setB)) - len(setB.difference(setA))

if __name__ == '__main__':
	#Setting up the argument parser
	parser = args.ArgumentParser(description="This script calculates all locally similar intervals inside sketches of two sequences.")
	parser.add_argument('-a', metavar='SeqA', type=str, required=True, help="The first sequence")
	parser.add_argument('-b', metavar='SeqB', type=str, required=True, help="The second sequence")
	parser.add_argument('-k', metavar='KmerLen', type=int, default=9, help="The k-mer length to use")
	parser.add_argument('-r', metavar='HashRatio', type=float, default=0.1, help="The ratio of hash values from the set of all possible hash " \
		"values to be included into a sketch")

	arguments = parser.parse_args()

	#Calculate hash threshold
	t = floor(((4 ** arguments.k) - 1) * arguments.r)
	#Calculate sketches
	skA = calcSketch(arguments.a, arguments.k, t)
	skB = calcSketch(arguments.b, arguments.k, t)

	#Calculate maximum similarities
	results = []

	##Attempt to prune the search space##

	#We start with the largest intervals possible
	intCandidates = [(0, 0, len(skA) - 1, 0, len(skB) - 1)]

	#Iterate over interval canditates as long as there are any
	while len(intCandidates) > 0:
		#Get next candidate
		candidate = intCandidates.pop(0)
		#Score the candidate
		candidate[0] = calcSim(skA, candidate[1], candidate[2], skB, candidate[3], candidate[4])

		#If score is positive try to add it to results
		if candidate[0] > 0:
			addRes(results, candidate)

		#Get interval sizes
		lIntSize = candidate[2] - candidate[1] + 1
		rIntSize = candidate[4] - candidate[3] + 1

		#Generate new candidates from old candidate except if we have found the maximum possible score already
		if candidate[0] < min(lIntSize, rIntSize) and max(lIntSize, rIntSize) > 1:
			if lIntSize == rIntSize:
				intCandidates.append((0, candidate[1], candidate[2], candidate[3] + 1, candidate[4]))
				intCandidates.append((0, candidate[1], candidate[2], candidate[3], candidate[4] - 1))
			else:
				intCandidates.append((0, candidate[1] + 1, candidate[2], candidate[3], candidate[4]))
				intCandidates.append((0, candidate[1], candidate[2] - 1, candidate[3], candidate[4]))

	##Most naive approach##

	#Iterate over all possible intervals in seqA
	for s in range(len(skA)):
		for e in range(len(skA)):
			if s <= e:
				#Iterate over all possible intervals in seqB
				for u in range(len(skB)):
					for f in range(len(skB)):
						if u <= f:
							score = calcSim(skA, s, e, skB, u, f)

							#We are only interested in positive results
							if score > 0:
								addRes = True

								#Check if there is already a result overlapping this one
								for r in results:
									if (s <= r[1] and r[1] <= e) or (u <= r[3] and r[3] <= f) or (r[1] <= s and s <= r[2]) or \
									(r[3] <= u and u <= r[4]):
										if score < r[0]:
											addRes = False
											break
										
										if score == r[0] and (e - s + 1) + (f - u + 1) <= (r[2] - r[1] + 1) + (r[4] - r[3] + 1):
											addRes = False
											break
										
										results.remove(r)

								if addRes:
									results.append((score, s, e, u, f))

	print(results)
