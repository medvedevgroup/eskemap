#!/usr/bin/env python3

import sys

#This script transforms the output of different implementations of the DP algorithm into a unified format.

#Parse file
fh = open(sys.argv[1], 'r'):
firstLine = fh.readline()

#Check which format we are dealing with
if firstLine.startswith("Pair"):
	for l in fh.readlines():
		l = l.strip()

		if l.startswith('['):
			fstSplit = l.split(' ')
			i, j = fstSplit[0].replace('[', "").replace("]:", "").split(',')

			print(f"{i} {j} {fstSplit[1]}")
else:
	for l in firstLine + fh.readlines():
		l = l.strip()
		elems = l.split(' ')

		print(f"{elems[1]} {elems[3]} {elems[5]}")
