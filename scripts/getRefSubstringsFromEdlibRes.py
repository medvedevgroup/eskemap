#!/usr/bin/env python3

import argparse as args
from Bio import SeqIO
from sys import stderr
from sys import maxsize
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#This script returns a FASTA file containing substrings of a given sequence which correspond to reported mapping positions within a 
#given merged result file of edlib

#Setting up the argument parser
parser = args.ArgumentParser(description="This script creates a FASTA file containing substrings of a given reference sequence " + \
	"which correspond to a mapping location as specified in a given merged edlib mapping file.")
parser.add_argument('-r', metavar='ref', type=args.FileType('r'), required=True, help="Reference sequence file in FASTA format")
parser.add_argument('-m', metavar='maps', type=args.FileType('r'), required=True, help="Mapping file")
parser.add_argument('-t', metavar='thres', type=int, required=True, help="Only consider reads with at most this many mappings")
parser.add_argument('-o', metavar='ofile', type=args.FileType('w'), required=True, help="Output file name")

arguments = parser.parse_args()
#Load reference sequence
refRecs = [r for r in SeqIO.parse(arguments.r, "fasta")]

#Drop a warning that in case of multiple reference sequences we only consider the first one
if len(refRecs) > 1:
	print("Warning: Multiple reference sequences found. Only first one is considered", file=stderr)

#Get the reference sequence of interest
refSeq = str(refRecs[0].seq)
#A list to save substrings to output
substringRecs = []
#A list of mapping coordinates
coords = []

#Read mapping file
for l in arguments.m:
	#Check if mappings for a new read start
	if l.endswith(".er\n"):
		#Count number of mappings for last read
		nbMappings = len(coords)

		#Check if mappings of previous read should be added to output
		if nbMappings > 0  and nbMappings <= arguments.t:
			#Iterate over mapping coordinates
			for c in coords:
				#Create sequence record for corresponding substring
				substringRecs.append(SeqRecord(Seq(refSeq[c[0]:c[1]+1]), id=f"s_{rid}:ref{c[0]}-{c[1]}"))

		#Clear coordinate list
		coords = []
		#Parse read id
		rid = l.split("_ri")[1].split(".er")[0]
	else:
		#Parse mapping coordinates
		coords.append([int(i) for i in l.split(' ') if i.isdecimal()][:2])

#Also do the check for the last mapping
if nbMappings > 0  and nbMappings <= arguments.t:
	#Iterate over mapping coordinates
	for c in coords:
		#Create sequence record for corresponding substring
		substringRecs.append(SeqRecord(Seq(refSeq[c[0]:c[1]+1]), id=f"s_{rid}:ref{c[0]}-{c[1]}"))

#Write data to file
SeqIO.write(substringRecs, arguments.o, "fasta")
