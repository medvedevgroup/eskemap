#!/usr/bin/env python3

import argparse as args
from math import inf, floor
from CalcLocInterSim import calcSketch
from sys import stdout
from os.path import exists
from Bio import SeqIO
from numpy import unique

#This script calculates the intersection similarity measure (as defined on June 16 2022 with a modification to treat duplications 
#from July 11 2022) between a given pattern and a text

if __name__ == '__main__':
	#Setting up the argument parser
	parser = args.ArgumentParser(description="This script calculates the intersection similarity between two sequences.")
	parser.add_argument('-p', metavar='Pattern', type=str, required=True, help="The pattern sequence (plain sequence or FASTA " + \
		"file)")
	parser.add_argument('-s', metavar='Text', type=str, required=True, help="The text sequence (plain sequence or FASTA file)")
	parser.add_argument('-k', metavar='KmerLen', type=int, default=9, help="The k-mer length to use")
	parser.add_argument('-r', metavar='HashRatio', type=float, default=0.1, help="The ratio of hash values from the set of all " + \
		"possible hash values to be included into a sketch")
	parser.add_argument('-c', metavar='CommonWeight', type=int, default=1, help="Weight to reward common hashes (default 1)")
	parser.add_argument('-u', metavar='UniqueWeight', type=int, default=1, help="Weight to punish unique hashes (default 1)")

	arguments = parser.parse_args()

	#Check if given weights are positive
	if arguments.c < 1 or arguments.u < 1:
		print("ERROR: Weights must be positive!", file=stdout)
		exit(-1)

	#Calculate hash threshold
	ht = floor(((4 ** arguments.k) - 1) * arguments.r)

	#Calculate sketches from sequences
	if exists(arguments.p):
		patternSeq = [str(r.seq) for r in SeqIO.parse(open(arguments.p, 'r'), "fasta")][0]
	else:
		patternSeq = arguments.p

	pattern = calcSketch(patternSeq, arguments.k, ht)

	if exists(arguments.s):
		textSeq = [str(r.seq) for r in SeqIO.parse(open(arguments.s, 'r'), "fasta")][0]
	else:
		textSeq = arguments.s

	text = calcSketch(textSeq, arguments.k, ht)

	#Initialize score
	if text[0] in pattern:
		score = arguments.c - (len(pattern) - 1) * arguments.u
	else:
		score = -arguments.u * (1 + len(pattern))

	for i in range(1, len(text)):
		if text[i] in pattern:
			score += arguments.c + arguments.u
		else:
			score -= arguments.u

	print(score)
