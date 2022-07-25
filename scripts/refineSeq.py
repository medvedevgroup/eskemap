#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO, Seq
from Bio.Alphabet import _verify_alphabet
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from random import choice

#Constants
IUPAC_NON_SPEC = "RYSWKMBDHVN.-"
#NOTE: Gap symbols '.' and '-' are treated as if they could be replaced by any nucleotid to have a placeholder in this case
IUPAC_MEANING = {'R': "AG", 'Y': "CT", 'S': "GC", 'W': "AT", 'K': "GT", 'M': "AC", 'B': "CGT", 'D': "AGT", 'H': "ACT", 'V': "ACG", 'N': "ACGT", '.': "ACGT", '-': "ACGT"}

#Set argument parser
parser = argparse.ArgumentParser(description="This script refines input nucleotide sequences containing IUPAC characters producing sequences only consisting of explicit nucleotides (A,C,G and T) and an auxiliary file for trace back.")
parser.add_argument("seqFiles", help="FASTA files containing IUPAC characters to be refined", nargs='+', metavar="IUPAC_Seq_File")
parser.add_argument("-l", "--seqLen", help="Length of sequences in trace back file", type=int, default=31)
parser.add_argument("-o", "--outFile", help="Trace back file name to save IUPAC information", default="IUPAC_charTrcBck.fasta")
args = parser.parse_args()

#Open trace back file
# Format is
# >[Source file name] [origSeq] [ReplacedSeqVarInOrigFile]
# [ReplacedSeqVarInTraceFile]
tFile = open(args.outFile, 'w')

#Iterate over input files
for f in args.seqFiles:
	#Get all records
	records = SeqIO.parse(f, "fasta")
	#Initialize refined output string
	refOut = ""
	#Read in as FASTA and iterate over sequence records
	for rec in records:
		#Iterate over non ACGT IUPAC characters
		for c in IUPAC_NON_SPEC:
			#Search for IUPAC character
			while not rec.seq.find(c) < 0:
				#Find position of IUPAC character
				pos = str(rec.seq).index(c)
				#Get start position of last trace back sequence in record's sequence
				lstPos = len(rec.seq) - args.seqLen

				#Make sure we get a trace back sequence that is long enough
				if pos >= lstPos:
					#Get the record's sequence end
					origSeq = str(rec.seq[lstPos:])
				else:
					#Get trace back sequence
					origSeq = str(rec.seq[pos:][:args.seqLen])

				#This will be the refined sequence to be inserted in the original file
				refSeq = ""

				#Go through original sequence and refine it
				for n in origSeq:
					#Check if this base has to be refined
					if n in IUPAC_MEANING:
						#Refine it
						refSeq += choice(IUPAC_MEANING[n])
					else:
						#Take over old base
						refSeq += n

				#Find variant position
				varPos = origSeq.index(c)
				#Generate sequence variant for aux file
				varSeq = refSeq[:varPos] + choice(IUPAC_MEANING[origSeq[varPos]].replace(refSeq[varPos], '')) + refSeq[varPos + 1:]
				#Write aux file entry
				tFile.write(">" + f + ":" + rec.description + " orig:" + origSeq + " ref:" + refSeq + "\n" + varSeq + "\n")
				#Replace unrefined sequence
				rec.seq = Seq.Seq(str(rec.seq).replace(origSeq, refSeq, 1), IUPACUnambiguousDNA)

		#Output refined sequence
		refOut += '>%s\n%s\n' %(rec.description,rec.seq)

	#Update input file
	unRefFile = open(f, 'w')
	unRefFile.write(refOut)
	unRefFile.close()

tFile.close()