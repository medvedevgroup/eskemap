#!/usr/bin/env python3

import argparse as args
from random import choice, randrange, seed
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from sys import maxsize

#This script generates a random DNA sequence of specified length and saves it to file in FASTA format

#The nucleotide alphabet
NUCL_ALPHABET = "ACGT"

#Setting up the argument parser
parser = args.ArgumentParser(description="This script generates a random DNA sequence and writes it to file in FASTA format.")
parser.add_argument('-l', metavar='Length', type=int, required=True, help="Length of generated sequence")
parser.add_argument('-o', metavar='Output', type=str, required=True, help="Name of output FASTA")
parser.add_argument('-s', metavar='Seed', type=int, default=randrange(maxsize), help="Random seed to use")

arguments = parser.parse_args()

#Initialize random seed
seed(arguments.s)

seq = ""

for i in range(arguments.l):
	seq += choice(NUCL_ALPHABET)

SeqIO.write([SeqRecord(Seq(seq), id="RandomSeq", description=f"Seed={arguments.s}")], open(arguments.o, 'w'), "fasta")
