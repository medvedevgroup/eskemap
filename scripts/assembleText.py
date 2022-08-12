#!/usr/bin/env python3

import argparse as args
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#This script assembles a text taking some sequence as base and injecting a given set of other sequences into it. Injected sequences 
#are spread evenly over the base sequence.

#Setting up the argument parser
parser = args.ArgumentParser(description="This script assembles a text from some base and a set of sequences to inject.")
parser.add_argument('-b', metavar='TextBase', type=str, required=True, help="File containing the text base")
parser.add_argument('-i', metavar='InsSeq', type=str, required=True, nargs='+', help="Files containing sequences to insert")
parser.add_argument('-o', metavar='Output', type=str, required=True, help="Name of output FASTA")

arguments = parser.parse_args()

#Load text base
textbase = [r.seq for r in SeqIO.parse(open(arguments.b, 'r'), "fasta")][0]

#Load sequences to inject
injSeqs = []

for f in arguments.i:
	injSeqs += [r.seq for r in SeqIO.parse(open(f, 'r'), "fasta")]

#Construct text
basePartLen = int(len(textbase) / (len(injSeqs) + 1))
text = textbase[:basePartLen]

for i in range(len(injSeqs)):
	text += injSeqs[i] + textbase[basePartLen + i * basePartLen: (i + 2) * basePartLen]

rid = arguments.b.split("_rid")[1].split(".fasta")[0]

#Write text to file
SeqIO.write(SeqRecord(text, id="Text", description=f"rid={rid}"), open(arguments.o, 'w'), "fasta")
