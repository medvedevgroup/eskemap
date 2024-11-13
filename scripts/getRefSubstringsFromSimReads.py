#!/usr/bin/env python3

import argparse as args
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from sys import stderr

#This script returns a FASTA file containing the reference substrings of a given reference as specified in the headers of the 
#given read file

#Setting up the argument parser
parser = args.ArgumentParser(description="This script creates a FASTA file containing the reference substrings of a given refer" + \
    "ence as specified in the headers of the given read file")
parser.add_argument('-r', metavar='ref', type=args.FileType('r'), required=True, help="The reference sequence")
parser.add_argument('-d', metavar='reads', type=args.FileType('r'), required=True, help="The reads file")
parser.add_argument('-p', metavar='padding', type=float, default=1.0, help="Relative amount of surrounding sequence to include")
parser.add_argument('-o', metavar='ofile', type=str, required=True, help="Name of output file")

arguments = parser.parse_args()
#Load reference sequence
refRecs = [r for r in SeqIO.parse(arguments.r, "fasta")]

#Drop a warning that in case of multiple reference sequences we only consider the first one
if len(refRecs) > 1:
	print("Warning: Found multiple sequences in reference file. Only the first one is considered", file=stderr)

#Get the reference sequence of interest
refSeq = str(refRecs[0].seq)
#Get reference sequence coordinates from read file
refCoordStrings = [r.description.split(' ')[2] for r in SeqIO.parse(arguments.d, "fasta")]
substringRecs = []
#Set exit state
exitState = 0

#Iterate over all reads
for s in refCoordStrings:
    #Check if we are dealing with the correct reference
    if refRecs[0].id != s.split(':')[0]:
        print("ERROR: Reference name from read header does not match the given reference's name!", file=stderr)
        exitState = -1
        break

    #Extract start and end coordinates
    start, end = [int(p) for p in s.split(':')[1].split('-')]
    #Calculate number of padding bases at both flanks
    nbPadBases = int((arguments.p * (end - start + 1)) / 2 + 0.5)
    #Adjust start and end depending on padding parameter
    start = max(0, start - nbPadBases)
    end = min(len(refSeq), end + nbPadBases)
    substringRecs.append(SeqRecord(Seq(refSeq[start: end + 1]), id=refRecs[0].id, description=f"substring: {start}-{end}"))

#Write sequences to disc
SeqIO.write(substringRecs, open(arguments.o, 'w'), "fasta")
exit(exitState)
