#!/usr/bin/env python3

import argparse as args
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from sys import stderr

#This script extracts the i-th entry from a multiple FASTA file

#Setting up the argument parser
parser = args.ArgumentParser(description="This script extracts the i-th entry (starting from 1) from a multiple FASTA file.")
parser.add_argument('-s', metavar='Source', type=str, required=True, help="File to extract from")
parser.add_argument('-i', metavar='EtrIdx', type=int, required=True, help="Index of entry to extract")
parser.add_argument('-o', metavar='Ofile', type=str, required=True, help="Output file name")

arguments = parser.parse_args()

entries = SeqIO.parse(open(arguments.s, 'r'), "fasta")

if arguments.i > len(entries):
	print("ERROR: i-th entry does not exist in input file!", file=stderr)
	exit(-1)

SeqIO.write(entries[arguments.i - 1], open(arguments.o, 'w'), "fasta")
