#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO

#Set argument parser
parser = argparse.ArgumentParser(description="This script takes a refined nucleotide sequence and a trace back file as input and undoes the refinement.")
parser.add_argument("traceFile", help="Trace back file containing applied refinement infos")
parser.add_argument("-s", "--source", help="FASTA files containing the refined input sequences", nargs='+')
parser.add_argument("-v", "--variants", help="File containing refined variants")
args = parser.parse_args()

#Check if input suffices to do anything
if not args.source and not args.variants:
	print("ERROR: Missing input files\nEither refined input sequences or called variants are required for undoing refinements\nCannot do anything", file=sys.stderr)
	exit(1)

#Dictionary to store traces: {(seqFile, recordID): [(origSeq, refSeq), ...], ...}
traces = {}

#Read trace back file
for trace in SeqIO.parse(args.traceFile, "fasta"):
	#Testing
	# print(trace.description)

	#Parse record id
	#Split location and sequence information
	loc, seq = trace.description.split(" orig:")
	#Split location infos
	fileName, seqHeader = loc.split(':')
	#Split sequence infos
	seqInfos = seq.split(" ref:")

	#Check if there exists an entry for this sequence already
	if (fileName, seqHeader) in traces:
		#Append sequence infos
		traces[(fileName, seqHeader)].append(seqInfos)
	else:
		#Make new entry
		traces[(fileName, seqHeader)] = [seqInfos]

#Testing
# print(traces)
# exit(0)

#Check if there are source files to unrefine
if args.source:
	#Iterate over source files
	for f in args.source:
		records = {}

		#Parse all records in file
		for r in SeqIO.parse(f, "fasta"):
			#Add record
			records[r.description] = str(r.seq)

			#Check if there are refinements for this record
			if (f,r.description) in traces:
				#Get refinements performed for this record
				refs = list(traces[(f,r.description)])

				#Iterate over refinements for undoing them starting with the last one applied
				while refs != []:
					#Get a refinement
					ref = refs.pop()

					#Testing
					# print(records[r.description].replace("A", "B"))
					# print(ref[0])
					# print(records[r.description])

					#Replace refined sequence part by original one
					records[r.description] = records[r.description].replace(ref[1], ref[0])

					#Testing
					# print(records[r.description])

		#Update file
		sFile = open(f, 'w')
		#Iterate over records
		for r in records:
			#Write record
			sFile.write(">%s\n%s\n" %(r, records[r]))

		#Close file
		sFile.close()

#Check if there is a variant file to unrefine
if args.variants:
	#Initialize output string
	outStr = ""
	#Read variant file
	vFile = open(args.variants, 'r')
	#Initialize dictionary to store variant sequences
	varSeqs = {}
	#Last core unitig sequence
	lstCore = ""

	#Iterate over lines
	for l in vFile:
		#Remove line break
		l = l.rstrip('\n')

		#Check if current line is a core sequence
		if l.startswith("Core unitig: "):
			#Extract sequence
			lstCore = l.split(' ')[2]
			#Get core sequence length
			lstCoreLen = len(lstCore)

			#Add core sequence to all variants
			for v in varSeqs:
				varSeqs[v][0] += lstCore
				varSeqs[v][1].append(lstCoreLen)
		else:
			#Parse line
			variant, seq = l.split(':')

			#Check if we have seen this variant sequence before
			if variant in varSeqs:
				#Add sequence
				varSeqs[variant][0] += seq
				varSeqs[variant][1].append(len(seq))
			#Check if we have seen the first core unitig already
			elif lstCore != "":
				#If a core unitig already exists save that this sequence variant has no sequence left from it
				varSeqs[variant] = ["", [0]]
				#Add last core sequence and variant sequence
				varSeqs[variant][0] += lstCore + seq
				varSeqs[variant][1].append(lstCoreLen)
				varSeqs[variant][1].append(len(seq))
			else:
				#Add new entry
				varSeqs[variant] = [seq, [len(seq)]]

	#Close file
	vFile.close()

	#Go through traces
	for t in traces:
		#Iterate over refinements
		while traces[t] != []:
			#Get last applied refinement
			ref = traces[t].pop()
			#Undo refinement
			varSeqs[t[0]][0] = varSeqs[t[0]][0].replace(ref[1], ref[0])

	#Open variant file again
	#vFile = open(args.variants, 'w')
	#A flag to indicate that we are currently outputting variant sequences in front of the first core unitig
	inFrntOfFstCore = True
	#A flag indicating what kind of sequence we are going to report next (core unitig or variant sequence)
	repVarSeq = True

	#Make sure that all sequences have been completely reported
	while max([len(v[1]) for v in varSeqs.values()]) > 0:
		#Iterate over variants
		for v in varSeqs:
			#Check if there is no sequence left for this variant to be reported
			if varSeqs[v][1] == []:
				continue

			#Get length of next substring to output
			substrLen = varSeqs[v][1][0]
			#Remove that length from list
			varSeqs[v][1].remove(substrLen)

			#Check if we are going to report sequence variants
			if repVarSeq:
				#If we have not seen the first core unitig yet and there is no variant sequence to report we output nothing
				if inFrntOfFstCore and substrLen == 0:
					continue
		
				#Report variant sequence
				print(v + ':' + varSeqs[v][0][:substrLen])
			else:
				#Get the core sequence to report
				coreSeq = varSeqs[v][0][:substrLen]

			#Remove substring from sequence
			varSeqs[v][0] = varSeqs[v][0][substrLen:]

		#Check if we still owe to report a core sequence
		if not repVarSeq:
			#Report core sequence
			print("Core unitig:", coreSeq)

		#We are going to report the opposite kind of sequence next
		repVarSeq = not repVarSeq
		#All variant sequences in front of first core are reported
		inFrntOfFstCore = False

	#Close file
	# vFile.close()
