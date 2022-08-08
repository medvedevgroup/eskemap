#!/usr/bin/env python3

import sys

#This script transforms the output of different implementations of the DP algorithm into a unified format.

#Parse file

#Check which format we are dealing with
for l in open(sys.argv[1], 'r'):
	l = l.strip()

	#Script implementation format
	if l.startswith('['):
		fstSplit = l.split(' ')
		i, j = fstSplit[0].replace('[', "").replace("]:", "").split(',')

		print(f"{i} {j} {fstSplit[1]}")

	#C++ implementation format
	if l.startswith("i: "):
		elems = l.split(' ')

		print(f"{elems[1]} {elems[3]} {elems[5]}")
