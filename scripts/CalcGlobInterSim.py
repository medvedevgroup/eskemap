#!/usr/bin/env python3

import argparse as args
from math import floor
from sys import stdout
from os.path import exists
from Bio import SeqIO
from numpy import unique

#This script calculates the intersection similarity measure (as defined on June 16 2022 with a modification to treat duplications 
#from July 11 2022) between a given pair of sketches

if __name__ == '__main__':
	#Setting up the argument parser
	parser = args.ArgumentParser(description="This script calculates the global intersection similarity between a pair of " + \
		"sketches.")
	parser.add_argument('-p', metavar='Pair', type=str, required=True, help="The sketch pair (.sk format) file)")
	parser.add_argument('-c', metavar='CommonWeight', type=int, default=1, help="Weight to reward common hashes (default 1)")
	parser.add_argument('-u', metavar='UniqueWeight', type=int, default=1, help="Weight to punish unique hashes (default 1)")

	arguments = parser.parse_args()

	#Check if given weights are positive
	if arguments.c < 1 or arguments.u < 1:
		print("ERROR: Weights must be positive!", file=stdout)
		exit(-1)

	#Load the sketches
	sketches = []

	for l in open(arguments.p, 'r'):
		if l.startswith('>'):
			continue

		sketches.append(l.split(' '))

	#Initialize score
	if sketches[0][0] in sketches[1]:
		score = arguments.c - (len(sketches[1]) - 1) * arguments.u
	else:
		score = -arguments.u * (1 + len(sketches[1]))

	for i in range(1, len(sketches[0])):
		if sketches[0][i] in sketches[1]:
			score += arguments.c + arguments.u
		else:
			score -= arguments.u

	print(score)
