#!/usr/bin/env python3

import argparse as args
from MutateSeq import mutateSeq
from random import randrange, randint, seed
from sys import maxsize
from math import ceil
from numpy import random
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

#This script simulates reads from a given reference sequence using a mutational model as described in the manuscript

if __name__ == '__main__':
	#Setting up the argument parser
	parser = args.ArgumentParser(description="This script simulates reads from a given reference sequence.")
	parser.add_argument('-dp', metavar='Depth', type=int, required=True, help="Sequencing depth to simulate")
	parser.add_argument('-lmn', metavar='LenMin', type=int, required=True, help="Minimum read length")
	parser.add_argument('-lmx', metavar='LenMax', type=int, required=True, help="Maximum read length")
	parser.add_argument('-lavg', metavar='LenMean', type=float, required=True, help="Average read length")
	parser.add_argument('-ls', metavar='LenStd', type=float, required=True, help="Standart deviation of read lengths")
	parser.add_argument('-sr', metavar='SubRate', type=float, required=True, help="Substitution rate")
	parser.add_argument('-dr', metavar='DelRate', type=float, required=True, help="Deletion rate")
	parser.add_argument('-ir', metavar='InsLen', type=float, required=True, help="Mean insertion length")
	parser.add_argument('-sd', metavar='Seed', type=int, default=randrange(maxsize), help="Random seed to use")
	parser.add_argument('-r', metavar='Ref', type=str, required=True, help="Reference FASTA")
	parser.add_argument('-o', metavar='OutFile', type=str, required=True, help="Name of output file")

	arguments = parser.parse_args()
	#We assume that the reference sequence is a single continous sequence
	refSeq = [r for r in SeqIO.parse(open(arguments.r, 'r'), "fasta")][0]
	seed(arguments.sd)
	random.seed(randrange(2**32 - 1))
	rds = []

	#For each read
	for i in range(ceil((arguments.dp * len(refSeq.seq) / arguments.lavg))):
		rLen = -1

		#Determine a length
		while rLen < arguments.lmn or rLen > arguments.lmx:
			rLen = int(random.gamma((arguments.lavg / arguments.ls)**2, scale=(arguments.ls**2) / arguments.lavg))

		#Determine start
		rs = randint(0, len(refSeq.seq) - rLen)
		#Get read sequence
		rdSeq = mutateSeq(str(refSeq.seq[rs:rs+rLen]), arguments.sr, arguments.dr, arguments.ir)
		rds.append(SeqRecord(Seq(rdSeq), id='s_' + str(i), description=f"ref: {refSeq.id}:{rs}-{rs+rLen-1} sr: {arguments.sr} ir: {arguments.ir}" + \
			f" dr: {arguments.dr} sd: {arguments.sd}"))

	#Output reads
	SeqIO.write(rds, open(arguments.o, 'w'), "fasta")
