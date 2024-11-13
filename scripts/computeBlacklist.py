#!/usr/bin/env python3

import argparse as args
from Bio.SeqIO import parse
from collections import deque
from Bio.Seq import Seq

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

#This script computes the minimizer sketches of all sequences in the given FASTA file and outputs the hashes of all k-mers from the 
# sketches which a occur more than o times, where o is a user specified parameter
if __name__ == '__main__':
    #Setting up the argument parser
    parser = args.ArgumentParser(description="This script computes the minimizer sketches of all given sequences and outputs th" + \
        "e hashes of k-mers occurring at least a certain number of times")
    parser.add_argument('-s', metavar="seqs", type=args.FileType('r'), required=True, help="The sequence file in FASTA format")
    parser.add_argument('-k', metavar="kLen", type=int, default=15, help="The k-mer length")
    parser.add_argument('-w', metavar="windowSize", type=int, default=10, help="Window size for minimizer calculation")
    parser.add_argument('-o', metavar="maxOccs", type=int, required=True, help="Occurrence threshold")
    arguments = parser.parse_args()
    kmerProfile = {}

    #Iterate over sequences 
    for r in parse(arguments.s, "fasta"):
        #Calculate minimizer sketch
        sketch = calcMiniSketch(r.seq, arguments.k, arguments.w)

        #Iterate over sketch and count k-mers
        for km in sketch:
            if not km[2] in kmerProfile:
                kmerProfile[km[2]] = 1
            else:
                kmerProfile[km[2]] += 1

    for km in kmerProfile:
        if kmerProfile[km] > arguments.o:
            print(km)
