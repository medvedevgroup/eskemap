#!/usr/bin/env python3

import argparse as args
from Bio import SeqIO
from sys import stderr
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#This script cuts out a specified string from a larger sequence

#Setting up the argument parser
parser = args.ArgumentParser(description="This script cuts out a specified string from a larger sequence.")
parser.add_argument('-i', metavar='Input', type=str, required=True, help="Sequence file to extract from")
parser.add_argument('-s', metavar='Start', type=int, required=True, help="Start of substring")
parser.add_argument('-e', metavar='End', type=int, required=True, help="End of substring")
parser.add_argument('-o', metavar='Ofile', type=str, required=True, help="Output file name")

arguments = parser.parse_args()

if arguments.s > arguments.e:
	print("ERROR: Invalid coordinates", file=stderr)

orig = [r for r in SeqIO.parse(open(arguments.i, 'r'), "fasta")][0]

if arguments.e >= len(orig.seq):
	print("ERROR: Input sequence is too short", file=stderr)

info = f"_substring{arguments.s}-{arguments.e}"
oSeq = SeqRecord(orig.seq[arguments.s:arguments.e + 1], id=orig.id + info, description=orig.description + info)
SeqIO.write(oSeq, open(arguments.o, 'w'), "fasta")
