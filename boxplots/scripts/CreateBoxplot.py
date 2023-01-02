#!/usr/bin/env python3

import argparse as args
from glob import glob
from matplotlib import pyplot as plt

#This script creates a boxplot using all available data for a given weight and ranges of mutation parameters. If no weight is given 
#Weighted Jaccard scores are used instead of linear scores.

#Setting up the argument parser
parser = args.ArgumentParser(description="This script creates a boxplot using all available data for a given weight and ranges " + \
	"of mutation parameters. If no weight is given Weighted Jaccard scores are used instead of linear scores.")
parser.add_argument('-l', metavar='SeqLen', type=str, required=True, help="Used sequence length")
parser.add_argument('-n', metavar='NumRep', type=str, required=True, help="Generated number of replications")
parser.add_argument('-w', metavar='Weight', type=str, help="Score weight used")
parser.add_argument('-s', metavar='SubRts', type=float, nargs='+', required=True, help="Substitution rates to include data for")
parser.add_argument('-d', metavar='DelRts', type=float, nargs='+', required=True, help="Deletion rates to include data for")
parser.add_argument('-i', metavar='InLens', type=float, nargs='+', required=True, help="Average insertion lengths to include " + \
	"data for")
parser.add_argument('-o', metavar='OutFl', type=str, default="out.pdf", help="Name of output PDF")

arguments = parser.parse_args()
scores = {}

for s in arguments.s:
	for d in arguments.d:
		for i in arguments.i:
			k = (s, d, i)
			scores[k] = []
			l = arguments.l
			pwc = f"scores/sketchPair_fmrandSeq_l{l}_ts*_smmutSeq_randSeq_l{l}_ts*_sr{s}_d{d}_i{i}_ms*_k15_r0.1"

			if arguments.w:
				w = arguments.w
				pwc += f"_w{w}.lscr"
			else:
				pwc += ".jscr"

			for f in glob(pwc):
				scores[k].append(float(open(f, 'r').readline()))

sortedKeys = sorted(scores.keys())
plt.boxplot([scores[k] for k in sortedKeys])

if arguments.w:
	plt.title("Linear score")
else:
	plt.title("Weighted Jaccard")

plt.xticks(range(1, len(sortedKeys) + 1), sortedKeys)
plt.xlabel("Mutation parameters")
plt.ylabel("Score")
plt.savefig(arguments.o, format="pdf")
