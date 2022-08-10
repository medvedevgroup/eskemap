#!/usr/bin/env python3

import argparse as args
from random import choice, randrange, seed, random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from sys import maxsize, stderr

#This script mutates a given DNA sequence in FASTA format according to a mutation process (described in the manuscript as of August 
#10 2022). Mutated sequence is outputted in FASTA format.

#Some constants#
#The nucleotide alphabet
NUCL_ALPHABET = "ACGT"
#Substitution alphabets
SUB_BASES = {'A': "CGT", 'C': "AGT", 'G': "ACT", 'T': "ACG"}

#Setting up the argument parser
parser = args.ArgumentParser(description="This script mutates a given DNA sequence and writes it to file in FASTA format.")
parser.add_argument('-m', metavar='SubRate', type=float, required=True, help="Substitution rate")
parser.add_argument('-d', metavar='DelRate', type=float, required=True, help="Deletion rate")
parser.add_argument('-i', metavar='InsLen', type=float, required=True, help="Mean insertion length")
parser.add_argument('-t', metavar='Template', type=str, required=True, help="File name of template sequence to mutate")
parser.add_argument('-o', metavar='Output', type=str, required=True, help="Name of output FASTA")
parser.add_argument('-s', metavar='Seed', type=int, default=randrange(maxsize), help="Random seed to use")

arguments = parser.parse_args()

#Check if rates and insertion length were chosen reasonable
if arguments.m < 0 or arguments.m > 1:
	print("ERROR: Substitution rate should be between 0 and 1", file=stderr)

if arguments.d < 0 or arguments.d > 1:
	print("ERROR: Deletion rate should be between 0 and 1", file=stderr)

if arguments.i < 0:
	print("ERROR: Mean insertion length should be positive", file=stderr)

#Load template
templateRecord = list(SeqIO.parse(open(arguments.t, 'r'), "fasta"))[0]
seq = [b for b in templateRecord.seq]
#Initialize random seed
seed(arguments.s)

#Apply substitutions and deletions
for i in range(len(seq)):
	#Throw the dice
	outcome = random()

	#Select an action
	if outcome <= arguments.m:
		seq[i] = choice(SUB_BASES[seq[i]])
	elif outcome <= arguments.m + arguments.d:
		#This means even if a base is deleted we will allow insertions to appear before and after its previous place
		seq[i] = ""

#Apply insertions
for i in range(1, len(seq) + 1):
	while random() > (1 / (1 + arguments.d)):
		seq.insert(i, choice(NUCL_ALPHABET))

#Write mutated sequence to file
SeqIO.write([SeqRecord(Seq("".join(seq)), id="MutatedSeq", description=f"{templateRecord.id}_{templateRecord.description}_" + \
	f"MutationSeed={arguments.s}")], open(arguments.o, 'w'), "fasta")
