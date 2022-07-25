#!/usr/bin/env python3

import argparse as args
from math import inf, floor
from CalcLocInterSim import calcSketch
from sys import stdout
from os.path import exists
from Bio import SeqIO

#This script outputs all t-homologies between two given sequences considering the intersection similarity measure (as defined on June 16 2022) and 
#using a dynamic programming algorithm 

#This function reports a t-homology
def reportThomology(intStart, intEnd, score, pattern, text):
	print(f"[{intStart},{intEnd}]:", score)
	# print("p:", pattern)
	# print("t[i,j]:", text[intStart:intEnd + 1])

if __name__ == '__main__':
	#Setting up the argument parser
	parser = args.ArgumentParser(description="This script searches for t-homologies between two sequences.")
	parser.add_argument('-p', metavar='Pattern', type=str, required=True, help="The pattern sequence (plain sequence or FASTA file)")
	parser.add_argument('-s', metavar='Text', type=str, required=True, help="The text sequence (plain sequence or FASTA file)")
	parser.add_argument('-t', metavar='Tthres', type=int, help="The t-homology threshold. Default is negative length of pattern " + \
		"sketch times weight for unique hashes")
	parser.add_argument('-k', metavar='KmerLen', type=int, default=9, help="The k-mer length to use")
	parser.add_argument('-r', metavar='HashRatio', type=float, default=0.1, help="The ratio of hash values from the set of all possible hash " \
		"values to be included into a sketch")
	parser.add_argument('-c', metavar='CommonWeight', type=int, default=1, help="Weight to reward common hashes (default 1)")
	parser.add_argument('-u', metavar='UniqueWeight', type=int, default=1, help="Weight to punish unique hashes (default 1)")

	arguments = parser.parse_args()

	#Check if given weights are positive
	if arguments.c < 1 or arguments.u < 1:
		print("ERROR: Weights must be positive!", file=stdout)
		exit(-1)

	#Calculate hash threshold
	ht = floor(((4 ** arguments.k) - 1) * arguments.r)
	
	#Calculate sketches from sequences
	if exists(arguments.p):
		patternSeq = [str(r.seq) for r in SeqIO.parse(open(arguments.p, 'r'), "fasta")][0]
	else:
		patternSeq = arguments.p

	pattern = calcSketch(patternSeq, arguments.k, ht)

	if exists(arguments.s):
		textSeq = [str(r.seq) for r in SeqIO.parse(open(arguments.s, 'r'), "fasta")][0]
	else:
		textSeq = arguments.s

	text = calcSketch(textSeq, arguments.k, ht)
	#Define the score matrix
	scores = []
	#Define occ_p map
	occp = {}
	#Define last position array in order to answer queries on u(s[i,j])
	lastPos = []

	#Set threshold correctly
	if not arguments.t:
		THRES = -arguments.u * len(pattern)
	else:
		THRES = arguments.t

	#Fill occ_p map
	for h in pattern:
		if h in occp:
			occp[h] += 1
		else:
			occp[h] = 1

	#Fill lastPos array
	seenHashes = {}

	for i in range(len(text)):
		if text[i] in seenHashes:
			lastPos.append(seenHashes[text[i]])
			seenHashes[text[i]] = i
		else:
			lastPos.append(-1)
			seenHashes[text[i]] = i

	#Testing
	# print("Pattern:",  pattern)
	# print("Text:", text)

	#Fill score matrix
	for j in range(len(text)):
		#Add a new list to store all scores of this column
		scores.append([])

		for i in range(len(text)):
			if j < i:
				continue
			elif i == j:
				#Calculate initial score
				if text[j] in occp:
					#Testing
					#print("Option 1/5")

					scores[j].append([i, arguments.c * (occp[text[j]] + 1) - arguments.u * (len(pattern) - occp[text[j]])])
				else:
					#Testing
					#print("Option 2/5")

					scores[j].append([i, arguments.u * (-1 - len(pattern))])
			else:
				#Update score
				scores[j].append(list(scores[j - 1][i]))

				if not text[j] in occp:
					#Testing
					#print("Option 3/5")

					scores[j][-1][1] -= arguments.u
				elif lastPos[j] >= i:
					#Testing
					#print("Option 4/5")

					scores[j][-1][1] += arguments.c
				else:
					#Testing
					#print("Option 5/5")

					scores[j][-1][1] += arguments.c * 2 + arguments.u * occp[text[j]]

	#Testing
	# print("scores:", scores)
	# exit(0)

	#Walk through score matrix and find maximal t-homologies
	maxScores = {}

	for j in range(len(text) - 1, -1, -1):
		for i in range(len(text)):
			#We only have values for the upper half of the matrix
			if i > j:
				continue

			#Check if we are in the first row
			if i == 0:
				#First and last hash in substring need to be shared
				if text[i] in occp and text[j] in occp:
					#Check if we have already seen a relevant maximum
					if i in maxScores:
						#Check if we have found a maximal t-homology
						if scores[j][i][1] > maxScores[i]:
							#Update max
							maxScores[i] = scores[j][i][1]
							#Output t-homology
							reportThomology(i, j, scores[j][i][1], pattern, text)
					else:
						#Check if we have found a t-homology
						if scores[j][i][1] >= THRES:
							#Save new max
							maxScores[i] = scores[j][i][1]
							#Output t-homology
							reportThomology(i, j, scores[j][i][1], pattern, text)
			else:
				if i - 1 in maxScores and i in maxScores:
					#Take over maximum
					maxScores[i] = max(maxScores[i - 1], maxScores[i])
				#Check if we have only found a maximum in a preceding row
				elif i - 1 in maxScores:
					#Take over maximum
					maxScores[i] = maxScores[i - 1]

				#Check if we have found a t-homology
				if  text[i] in occp and text[j] in occp:
					if i in maxScores:
						if scores[j][i][1] > maxScores[i]:
							#Update max
							maxScores[i] = scores[j][i][1]
							#Output t-homology
							reportThomology(i, j, scores[j][i][1], pattern, text)
					else:
						if scores[j][i][1] >= THRES:
							#Update max
							maxScores[i] = scores[j][i][1]
							#Output t-homology
							reportThomology(i, j, scores[j][i][1], pattern, text)

			''' This is (hopefully) only a more complicated version of the above
			#Check if we are dealing with the last column
			if j == len(text) - 1:
				#The upper, right corner is an exceptional case
				if i == 0:
					#Check if we are dealing with a t-homology
					if scores[i][j] >= THRES:
						#Save new max
						maxScores[i] = scores[i][j]
						#Output t-homology
						reportThomology(i, j, scores[i][j], pattern, text)
				else:
					#Check if we have already found a maximum
					if i - 1 in maxScores:
						#In the last column a max for the last row is also a max for the current row
						maxScores[i] = maxScores[i - 1]
				
						#Check if we are dealing with a maximum t-homology
						if scores[i][j] > maxScores[i]:
							#Save new max
							maxScores[i] = scores[i][j]
							#Output t-homology
							reportThomology(i, j, scores[i][j], pattern, text)
					#Check if we have found the first t-homology
					elif scores[i][j] >= THRES:
						maxScores[i] = scores[i][j]
						#Output t-homology
						reportThomology(i, j, scores[i][j], pattern, text)
			else:
				#Check if we are dealing with the first row
				if i == 0:
					#Check if we have already seen a relevant maximum
					if i in maxScores:
						#Check if we have found a maximal t-homology
						if scores[i][j] > maxScores[i]:
							#Update max
							maxScores[i] = scores[i][j]
							#Output t-homology
							reportThomology(i, j, scores[i][j], pattern, text)
					else:
						#Check if we have found a t-homology
						if scores[i][j] >= THRES:
							#Save new max
							maxScores[i] = scores[i][j]
							#Output t-homology
							reportThomology(i, j, scores[i][j], pattern, text)
				#Check if we have already found a maximum in preceding rows and the last column
				else:
					if i - 1 in maxScores and i in maxScores:
						#Take over maximum
						maxScores[i] = max(maxScores[i - 1], maxScores[i])
					#Check if we have only found a maximum in a preceding row
					elif i - 1 in maxScores:
						#Take over maximum
						maxScores[i] = maxScores[i - 1]

					#Check if we have found a t-homology
					if i in maxScores and scores[i][j] > maxScores[i]:
						#Update max
						maxScores[i] = scores[i][j]
						#Output t-homology
						reportThomology(i, j, scores[i][j], pattern, text)
					elif scores[i][j] >= THRES:
						#Update max
						maxScores[i] = scores[i][j]
						#Output t-homology
						reportThomology(i, j, scores[i][j], pattern, text)
		'''
