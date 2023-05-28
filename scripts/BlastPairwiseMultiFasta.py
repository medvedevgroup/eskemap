#!/usr/bin/env python3

import argparse as args
from Bio import SeqIO
from os.path import dirname
from re import search
import subprocess
from sys import stderr
from os import remove
import shlex

parser = args.ArgumentParser(description="This script parses a multiple fasta file and runs pairwise blast for each record.")
parser.add_argument('-f', metavar='Fasta', type=str, required=True, help="Input FASTA")
parser.add_argument('-e', metavar='Evalue', type=float, default=10, help="E-value threshold to be used for BLAST")
parser.add_argument('-o', metavar='OutFile', type=str, required=True, help="Output file containing BLAST results")

arguments = parser.parse_args()
ofile = open(arguments.o, 'w')

for r in SeqIO.parse(open(arguments.f, 'r'), "fasta"):
	#Find out read id
	readid = search("s_[0-9]+", r.id + r.description)
	#Find out reference range
	rrange = search("ref[0-9]+-[0-9]+", r.id + r.description)

	#Testing
	# print(r.id)
	# if not readid:
	# 	print("Read id could not be found")
	# if not rrange:
	# 	print("Range could not be found")

	#Write single fasta file
	if dirname(arguments.f) == "":
		sfname = f"sub_{readid.group(0)}_{rrange.group(0)}.fasta"
	else:
		sfname = dirname(arguments.f) + f"/sub_{readid.group(0)}_{rrange.group(0)}.fasta"

	SeqIO.write(r, open(sfname, 'w'), "fasta")
	#Call BLAST
	# cmd = ['blastn', '-query', '../simulations/reads/t2thumanChrY_sr0.00010909090909090909_dr0.0009818181818181818_i0.000909090' + \
	# f'9090909091_sd7361077429744071834_lmn100_lmx1000000_lavg9000_ls7000_dp10_ri{readid.group(0)}.fasta', '-task', 'blastn', '-' + \
	# 'subject', sfname, '-outfmt', '"6 qacc qstart qend sacc sstart send evalue length pident nident mismatch positive gaps sstr' + \
	# 'and qcovhsp"']
	cmd = shlex.split("blastn -query ../simulations/reads/t2thumanChrY_sr0.00010909090909090909_dr0.0009818181818181818_i0.000909" + \
		f"0909090909091_sd7361077429744071834_lmn100_lmx1000000_lavg9000_ls7000_dp10_ri{readid.group(0).split('_')[1]}.fasta -tas" + \
		f"k blastn -subject {sfname} -outfmt '6 qacc qstart qend sacc sstart send evalue length pident nident mismatch positive g" + \
		f"aps sstrand qcovhsp' -evalue {arguments.e}")

	#Testing
	# print(cmd)

	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # call programm with defined out and err
	out, err = p.communicate() # open output channels

	#Testing
	# print(out.decode("utf-8"))

	if p.returncode != 0:
		print("ERROR: BLAST produced an error message:", err.decode("utf-8"), file=stderr)
		exit(-1)

	ofile.write(f"Results for: {readid.group(0)} {rrange.group(0)}\n")
	ofile.write(out.decode("utf-8"))
	#Remove single fasta file
	remove(sfname)

ofile.close()
