#!/usr/bin/env python3

import argparse as args
from Bio import SeqIO

#This script takes 

#Setting up the argument parser
parser = args.ArgumentParser(description="This script filters a set of reads based on the number of mapping results found by ed" + \
	"lib.")
parser.add_argument('-e', metavar='EdlibRes', type=args.FileType('r'), required=True, nargs='+', help="Edlib result file")
parser.add_argument('-r', metavar='Reads', type=args.FileType('r'), required=True, help="File containing the reads to filter")
parser.add_argument('-m', metavar='MaxRes', type=int, required=True, help="Maximum number of edlib results to filter")
parser.add_argument('-o', metavar='Output', type=args.FileType('w'), required=True, help="Name of output FASTA")

arguments = parser.parse_args()

#Load Edlib results
allEdlibRes = {}

for f in arguments.e:
    readId = int(f.name.split("_ri")[1].split(".er")[0])
    allEdlibRes[readId] = []
    
    for l in f:
        allEdlibRes[readId].append(l.strip().split(' '))

#Get a subset of all edlib results that is based on a maximum number of results per read
subSet = {}

for r in allEdlibRes:
	if len(allEdlibRes[r]) <= arguments.m:
		subSet[r] = list(allEdlibRes[r])

#Load read sequences
rdRecsDict = {rec.id: rec for rec in [rd for rd in SeqIO.parse(arguments.r, "fasta")]}

#Create a read file containing only those reads for which not too many edlib results have been found
SeqIO.write([rdRecsDict[f"s_{r}"] for r in subSet], arguments.o, "fasta")
