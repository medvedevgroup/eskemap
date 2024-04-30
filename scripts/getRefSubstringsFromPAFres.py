#!/usr/bin/env python3

import argparse as args
from Bio import SeqIO
from sys import stderr
import gzip
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#This function parses mapping coordinates from a mapping result file in compressed PAF format and returns it as a dictionary where 
#the keys are names of the mapped sequences and values are lists of pairs representing mapping positions with start and end coor-
#dinates
def parseCompressedPAF(filepath):
    resPerRead = {}
    
    for l in gzip.open(filepath, 'rt'):
        elems = l.strip().split('\t')
        
        # rdid = int(elems[0].split('_')[1])

        rdid = elems[0].split('_')[1]
        
        if rdid in resPerRead:
            resPerRead[rdid].append((int(elems[7]), int(elems[8]) - 1))
        else:
            resPerRead[rdid] = [(int(elems[7]), int(elems[8]) - 1)]
            
    return resPerRead

#This script returns a FASTA file containing substrings of a given sequence which correspond to reported mapping positions within a 
#given result file of some read mapping tool in compressed PAF format.
if __name__ == '__main__':
	#Setting up the argument parser
	parser = args.ArgumentParser(description="This script creates a FASTA file containing substrings of a given reference sequence " + \
		"which correspond to a mapping location as specified in a given sequence mapping file.")
	parser.add_argument('-r', metavar='ref', type=str, required=True, help="Reference sequence file in FASTA format")
	parser.add_argument('-p', metavar='maps', type=str, required=True, help="Mapping file (gzipped PAF format)")
	parser.add_argument('-f', metavar='ids', type=str, nargs='+', help="Names of mapped sequences to be considered")
	parser.add_argument('-o', metavar='ofile', type=str, required=True, help="Output file name")

	arguments = parser.parse_args()

	#Load reference sequence
	refRecs = [r for r in SeqIO.parse(open(arguments.r, 'r'), "fasta")]

	#Drop a warning that in case of multiple reference sequences we only consider the first one
	if len(refRecs) > 1:
		print("Warning: Multiple reference sequences found. Only first one is considered", file=stderr)

	#Get the reference sequence of interest
	refSeq = str(refRecs[0].seq)
	#Parse mapping file
	mappingCoordsPerRead = parseCompressedPAF(arguments.p)

	#If we have a specified list of sequence names we only need to output mapping positions for them
	if arguments.f:
		seqnames = arguments.f
	else:
		seqnames = mappingCoordsPerRead.keys()

	#A list to store sequence records for all substrings to be outputted
	substringRecs = []

	#Iterate over mapped sequences' names
	for n in seqnames:
		#Iterate over mapping coordinates
		for c in mappingCoordsPerRead[n]:
			substringRecs.append(SeqRecord(Seq(refSeq[c[0]:c[1]+1]), id=f"s_{n}:ref{c[0]}-{c[1]}"))

	#Write data to file
	SeqIO.write(substringRecs, open(arguments.o, 'w'), "fasta")
