#!/usr/bin/env python3

import argparse as args
from sys import stderr

#This script calculates the Weighted Jaccard score between a given 
#pair of sketches

if __name__ == '__main__':
	#Setting up the argument parser
	parser = args.ArgumentParser(description="This script calculates the Weighted Jaccard score between a pair of " + \
		"sketches.")
	parser.add_argument('-p', metavar='Pair', type=str, required=True, help="The sketch pair (.sk format) file)")

	arguments = parser.parse_args()

	#Read sketches and count k-mer occurrences
	occ = {}
	sketchCount = -1

	for l in open(arguments.p, 'r'):
		if l.startswith('>'):
			sketchCount += 1

			if sketchCount > 1:
				break

			continue

		l = l.strip()

		for k in l.split(' '):
			k = int(k)

			if not k in occ:
				occ[k] = [0, 0]

			occ[k][sketchCount] += 1.

	mins = 0
	maxs = 0

	for k in occ:
		mins += min(occ[k])
		maxs += max(occ[k])

	print(mins / maxs)
