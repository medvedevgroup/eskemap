#!/usr/bin/env python3

import argparse as args
from CalcLocInterSim import calcSketch
from MutateSeq import mutateSeq, NUCL_ALPHABET
from random import randrange, choice, seed
from sys import maxsize
from math import floor
from numpy import unique

#This script randomly generates a pair of sequences where one has been mutated according to a mutation process from the other. The 
#sequences have the special property that given a k-mer length and a min hash ratio their sketches contain no duplicates.

if __name__ == '__main__':
	#Setting up the argument parser
	parser = args.ArgumentParser(description="This script generates sketches of a sequence and its mutated copy that do not " + \
		"contain duplicates in .sk format.")
	parser.add_argument('-l', metavar='SeqLen', type=int, required=True, help="Length of sequences")
	parser.add_argument('-m', metavar='SubRate', type=float, required=True, help="Substitution rate")
	parser.add_argument('-d', metavar='DelRate', type=float, required=True, help="Deletion rate")
	parser.add_argument('-i', metavar='InsLen', type=float, required=True, help="Mean insertion length")
	parser.add_argument('-s', metavar='Seed', type=int, default=randrange(maxsize), help="Random seed to use")
	parser.add_argument('-k', metavar='KmerLen', type=int, default=9, help="The k-mer length to use")
	parser.add_argument('-r', metavar='HashRatio', type=float, default=0.1, help="The ratio of hash values from the set of all " + \
		"possible hash values to be included into a sketch")
	parser.add_argument('-n', action="store_true", help="Generate sketches that contain no duplicates")

	arguments = parser.parse_args()

	#Check if rates and insertion length were chosen reasonable
	if arguments.m < 0 or arguments.m > 1:
		print("ERROR: Substitution rate should be between 0 and 1", file=stderr)

	if arguments.d < 0 or arguments.d > 1:
		print("ERROR: Deletion rate should be between 0 and 1", file=stderr)

	if arguments.i < 0:
		print("ERROR: Mean insertion length should be positive", file=stderr)

	#Generate sequence whose sketch has no duplicates
	hasDuplicates = True
	seed(arguments.s)
	minHashThres = floor(((4 ** arguments.k) - 1) * arguments.r)

	while hasDuplicates:
		seq = ""

		for i in range(arguments.l):
			seq += choice(NUCL_ALPHABET)

		sketch = calcSketch(seq, arguments.k, minHashThres)

		if not arguments.n or len(sketch) == len(unique(sketch)):
			hasDuplicates = False

	#Mutate sequence and make sure its sketch still has no duplicates
	hasDuplicates = True

	while hasDuplicates:
		mutSeq = mutateSeq(seq, arguments.m, arguments.d, arguments.i)
		mutSeqSketch = calcSketch(mutSeq, arguments.k, minHashThres)

		if not arguments.n or len(mutSeqSketch) == len(unique(mutSeqSketch)):
			hasDuplicates = False

	#Output sketches
	if arguments.n:
		print(f">NoDuplicateSketchPair_{arguments.s} original template")
		print(' '.join([str(h) for h in sketch]))
		print(f">NoDuplicateSketchPair_{arguments.s} mutated template")
		print(' '.join([str(h) for h in mutSeqSketch]))
	else:
		print(f">SketchPair_{arguments.s} original template")
		print(' '.join([str(h) for h in sketch]))
		print(f">SketchPair_{arguments.s} mutated template")
		print(' '.join([str(h) for h in mutSeqSketch]))
