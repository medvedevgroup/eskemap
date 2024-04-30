#!/usr/bin/env python3

import argparse as args
from Bio import SeqIO
from sys import stderr
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import deque

NT_IN_BITS = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

#This function calculates the hash of a bitwise k-mer representation. The function is influenced by the code of "The
#minimizer Jaccard estimator is biased and inconsistent." from Belbasi et al. (function "minimap2_hash(seed,v,mask)"
#in file "minimap2_hash_uncompiled.py").
def getHash(kmer, mask):
    u = kmer & mask
    u = ((~u) + (u << 21)) & mask # u = (u<<21)-(u+1) = 77594587*u-1
    u = u ^ (u >> 24)
    u = ((u + (u << 3)) + (u << 8)) & mask # u *= 265
    u = u ^ (u >> 14)
    u = ((u + (u << 2)) + (u << 4)) & mask # u *= 21
    u = u ^ (u >> 28)
    u = (u + (u << 31)) & mask # u *= 2147483649

    return u

#This function calculates the minimizer sketch of a sequence. It is influenced by the code of "The minimizer Jaccard
#estimator is biased and inconsistent." from Belbasi et al. (function "winnowed_minimizers_linear(perm,windowSize)" 
#in file "winnowed_minimizers.py").
def calcMiniSketch(seq, k, w):
    sketch = []
    #A deque to store k-mer hashes inside the current window
    windowKmers = deque()
    mask = (4 ** k) - 1
    lastIdx = -1

    for i in range(len(seq) - k + 1):
        kmerBits = 0
        kmerBitsRevComp = 0
        windowBorder = i - (w - 1)
        
        #Get bit representation of k-mer
        for c in seq[i:i+k]:
            kmerBits = (kmerBits << 2) + NT_IN_BITS[c]

        #Get bit representation of k-mer's reverse complement
        for c in str(Seq(seq[i:i+k]).reverse_complement()):
            kmerBitsRevComp = (kmerBitsRevComp << 2) + NT_IN_BITS[c]

        #If a k-mer is its own reverse complement we skip it
        if kmerBits == kmerBitsRevComp:
            continue

        #Depending on which hash is smaller we consider either a k-mer or its reverse complement per position
        if kmerBits < kmerBitsRevComp:
            kmer = (i, kmerBits, getHash(kmerBits, mask))
        else:
            #A k-mer is a pair of k-mer's start position and its hash
            kmer = (i, kmerBitsRevComp, getHash(kmerBitsRevComp, mask))
            
        #Remove all k-mers with a hash value larger than the newly calculated one
        while (len(windowKmers) > 0) and (windowKmers[-1][2] > kmer[2]):
            windowKmers.pop()

        #Save new k-mer as window k-mer
        windowKmers.append(kmer)

        #Remove k-mer if it is not any longer inside the window
        while (len(windowKmers) > 0) and (windowKmers[0][0] < windowBorder):
            windowKmers.popleft()

        #As soon as we have seen a first full window of k-mers choose a minimizer
        if (windowBorder >= 0) and (len(windowKmers) > 0):      
            #We do not choose the same minimizer for a second time
            if lastIdx != windowKmers[0][0]:
                sketch.append((windowKmers[0][0]+k-1, windowKmers[0][1], windowKmers[0][2]))
                lastIdx = windowKmers[0][0]
                
            while len(windowKmers) > 1 and windowKmers[0][1] == windowKmers[1][1]:
                windowKmers.popleft()
                sketch.append((windowKmers[0][0]+k-1, windowKmers[0][1], windowKmers[0][2]))    
                lastIdx = windowKmers[0][0] 

    #If our sequence was too small to get a full window of k-mers to consider take the smallest one found so far
    if windowBorder < 0 and len(windowKmers) > 0:
        sketch.append((windowKmers[0][0]+k-1, windowKmers[0][1], windowKmers[0][2]))
        
        while len(windowKmers) > 1 and windowKmers[0][1] == windowKmers[1][1]:
            windowKmers.popleft()
            sketch.append((windowKmers[0][0]+k-1, windowKmers[0][1], windowKmers[0][2]))

    return sketch

#This script returns a FASTA file containing substrings of a given sequence which correspond to reported mapping positions within a 
#given result file of ESKEMAP
if __name__ == '__main__':
	#Setting up the argument parser
	parser = args.ArgumentParser(description="This script creates a FASTA file containing substrings of a given reference seque" + \
		"nce which correspond to a mapping location as specified in a given ESKEMAP mapping file.")
	parser.add_argument('-r', metavar='ref', type=str, required=True, help="Reference sequence file in FASTA format")
	parser.add_argument('-k', metavar='kLen', type=int, required=True, help="K-mer length used to calculate sketches for mapping")
	parser.add_argument('-w', metavar='wSize', type=int, required=True, help="Window size used to calculate sketches for mapping")
	parser.add_argument('-b', metavar='blist', type=args.FileType('r'), required=True, help="K-mer blacklist used to calcul" + \
		"ate sketches for mapping")
	parser.add_argument('-e', metavar='maps', type=str, required=True, help="Mapping file")
	parser.add_argument('-o', metavar='ofile', type=str, required=True, help="Output file name")

	arguments = parser.parse_args()
	#Load reference sequence
	refRecs = [r for r in SeqIO.parse(open(arguments.r, 'r'), "fasta")]

	#Drop a warning that in case of multiple reference sequences we only consider the first one
	if len(refRecs) > 1:
		print("Warning: Multiple reference sequences found. Only first one is considered", file=stderr)

	#Get the reference sequence of interest
	refSeq = str(refRecs[0].seq)
	#Load k-mer blacklist
	kmerBlackList = {}

	for l in arguments.b:
		kmerBlackList[int(l)] = None

	#Calculate reference sequence sketch
	refSketch = [k for k in calcMiniSketch(refSeq, arguments.k, arguments.w) if not k[2] in kmerBlackList]
	#Parse mapping file
	resPerRead = {}
	    
	for l in open(arguments.e, 'r'):
		l = l.strip()

		if l.startswith('s'):
			lastRdId = int(l.split('_')[1])
			resPerRead[lastRdId] = []
		else:
			cols = l.split(' ')
			s = refSketch[int(cols[1])][0] + 1 - arguments.k
			e = refSketch[int(cols[3])][0]
			resPerRead[lastRdId].append((s, e))

	#A list to store sequence records for all substrings to be outputted
	substringRecs = []

	#Iterate over mapping coordinates
	for r in resPerRead:
		for c in resPerRead[r]:
			substringRecs.append(SeqRecord(Seq(refSeq[c[0]:c[1]+1]), id=f"s_{r}:ref{c[0]}-{c[1]}"))

	#Write data to file
	SeqIO.write(substringRecs, open(arguments.o, 'w'), "fasta")
