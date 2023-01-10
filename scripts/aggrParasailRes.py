#!/usr/bin/env python3

import argparse as args

#This script takes multiple result files of parasail and aggregates them. Duplicates are removed.

#Setting up the argument parser
parser = args.ArgumentParser(description="This script takes multiple result files of parasail and aggregates them.")
parser.add_argument('-i', metavar='Results', type=str, required=True, nargs='+', help="Parasail result files")
parser.add_argument('-r', metavar='Recalc', type=bool, default=False, help="Recalculate borders according to cigar string")

arguments = parser.parse_args()

allRes = {}

for f in arguments.i:
	offset = int(f.split("_ra")[1].split('_')[0].split('-')[0])

	for l in open(f, 'r'):
		parts = l.strip().split(' ')

		for i in range(len(parts[:2])):
			parts[i] = int(parts[i]) + offset

		if parts[0] in allRes:
			if not parts[1] in allRes[parts[0]]:
				allRes[parts[0]][parts[1]] = parts[2]
			else:
				#Sanity check
				if parts[2] != allRes[parts[0]][parts[1]]:
					print("Warning: Equal results are not the same!")
		else:
			allRes[parts[0]] = {parts[1]: parts[2]}

for s in allRes:
	for e in allRes[s]:
		print(s, e, allRes[s][e])
