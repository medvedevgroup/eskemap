configfile: 'config.yaml'

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from os.path import exists
from random import randrange
from sys import maxsize
from glob import glob

NB_RANDSEQS = (config['nbCpys'] + 1) * config['nbSimSeqs']

def matchProgCallToSeed(wcs):
	samFiles = []
	sdToCall = {}

	for g in config['genomes']:
		for e in config['errorPatterns']:
			readSimLogFile = f"../simulations/reads/{g}_cR103_ep{e}.log"

			if exists(readSimLogFile):
				for l in open(readSimLogFile, 'r'):
					if l.startswith("seed"):
						seed = l.strip().split(' ')[2]
						break
			else:
				seed = str(randrange(maxsize))

			sdToCall[g, e] = seed

	for p in config['prog']:
		for g in config['genomes']:
			for e in config['errorPatterns']:
				samFiles.append(f"../simulations/{p}Res/{g}_cR103_ep{e}_s{sdToCall[g, e]}.sam.gz")

	return samFiles

def listUnifiedRes(wcs):
	return glob(f"../simulations/homologies/homologies_gn*_n*_p{wcs.prog}.uni")

def genAggList(wcs):
	return expand("../simulations/homologies/homologies_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{i}_d{d}_c{c}_u{u}_t{t}_n{n}_" + \
		"pCpp.txt",gn=wcs.gn, rn=wcs.rn, gl=wcs.gl, rl=wcs.rl, o=wcs.o, m=wcs.m, i=wcs.i, d=wcs.d, c=wcs.c, u=wcs.u, t=wcs.t, n=\
		range(int(config['nbSimSeqs'])))

rule all:
	input:
		expand("../simulations/homologies/homologies_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{m}_d{m}_{sp}_t0_pPy.txt", gn=\
			config['nbSimSeqs'], rn=NB_RANDSEQS, gl=config['geneLen'], rl=config['randSeqLen'], o=config['nbCpys'], m=\
			config['mutationRates'], sp=config['scoringPatterns']),
		expand("../simulations/homologies/homologies_gn{gn}_rn200_gl{gl}_rl1000_o1_m{m}_i{m}_d{m}_{sp}_t0_pPy.txt", gn=\
			config['nbSimSeqs'], gl=config['geneLen'], m=config['mutationRates'], sp=config['scoringPatterns']),
		expand("../simulations/nucmer/nucmerAlignments_gn{gn}_rn200_gl{gl}_rl1000_o1_m{m}_i{m}_d{m}_p{i}.coords", \
			gn=config['nbSimSeqs'], gl=config['geneLen'], m=config['mutationRates'], i=range(config['nbSimSeqs'])),
		expand("../simulations/nucmer/nucmerAlignments_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{m}_d{m}_p{i}.coords", gn=\
			config['nbSimSeqs'], rn=NB_RANDSEQS, gl=config['geneLen'], rl=config['randSeqLen'], o=config['nbCpys'], m=\
			config['mutationRates'], i=range(config['nbSimSeqs'])),
		expand("../simulations/nucmer/nucmerAlignments_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{m}_d{m}_maxmatch_p{i}.coords", gn=\
			config['nbSimSeqs'], rn=NB_RANDSEQS, gl=config['geneLen'], rl=config['randSeqLen'], o=config['nbCpys'], m=\
			config['mutationRates'], i=range(config['nbSimSeqs'])),
		expand("../simulations/homologies/homologies_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{m}_d{m}_{sp}_t0_n{n}_p{p}.uni", gn=\
			config['nbSimSeqs'], rn=NB_RANDSEQS, gl=config['geneLen'], rl=config['randSeqLen'], o=config['nbCpys'], m=\
			config['mutationRates'], sp=config['scoringPatterns'], n=range(100), p=config['implType']),
		expand("../simulations/homologies/homologies_gn{gn}_rn200_gl{gl}_rl1000_o1_m{m}_i{m}_d{m}_{sp}_t0_n{n}_p{p}.uni", gn=\
			config['nbSimSeqs'], gl=config['geneLen'], m=config['mutationRates'], sp=config['scoringPatterns'], n=range(100), p=\
			config['implType']),
		# expand("../simulations/allRes_p{p}.md5", p=config['implType'])
		# matchProgCallToSeed
		
# rule calcMd5:
# 	input:
# 		listUnifiedRes
# 	output:
# 		"../simulations/allRes_p{prog}.md5"
# 	shell:
# 		"md5sum {input} > {output}"

rule unifyRes:
	input:
		"../simulations/homologies/homologies_gn{desc}.txt"
	output:
		"../simulations/homologies/homologies_gn{desc}.uni"
	priority:
		1
	shell:
		"python3 scripts/unifyRes.py {input} > {output}"

rule buildBWAindex:
	input:
		"../simulations/genomes/{genome}.fasta"
	output:
		temp(expand("../simulations/genomes/{g}.fasta.{s}", g="{genome}", s=config['bwaIdxFileSufs']))
	shell:
		"bwa index {input}"

rule runBWAmem:
	input:
		ref = "../simulations/genomes/{genome}.fasta",
		idx = expand("../simulations/genomes/{g}.fasta.{s}", g="{genome}", s=config['bwaIdxFileSufs']),
		rds = "../simulations/reads/{genome}_c{chem}_ep{dRat}_s{seed}.fastq.gz"
	output:
		"../simulations/bwamemRes/{genome}_c{chem}_ep{dRat}_s{seed}.sam.gz"
	shell:
		"bwa mem {input.ref} {input.rds} | gzip -3 > {output}"

rule runMinimap2:
	input:
		ref = "../simulations/genomes/{genome}.fasta",
		qry = "../simulations/reads/{genome}_c{chem}_ep{dRat}_s{seed}.fastq.gz"
	output:
		"../simulations/minimap2Res/{genome}_c{chem}_ep{dRat}_s{seed}.sam.gz"
	shell:
		"minimap2 -ax map-ont {input.ref} {input.qry} | gzip -3 > {output}"

rule reproduceSimulation:
	input:
		ref = "../simulations/genomes/{genome}.fasta",
		model = "../software/pbsim2/data/{chem}.model"
	params:
		dr = "{dRat}",
		sd = "{seed}"
	output:
		rds = temp("../simulations/reads/{genome}_c{chem}_ep{dRat}_s{seed}.fastq.gz"),
		maf = temp("../simulations/reads/{genome}_c{chem}_ep{dRat}_s{seed}.maf.gz"),
		log = "../simulations/reads/{genome}_c{chem}_ep{dRat}_s{seed}.log"
	shell:
		"pbsim --hmm_model {input.model} --difference-ratio {params.dr} --prefix $(echo {output.log} | sed 's/.log//g') --depth" + \
		" 10.0 --seed {params.sd} {input.ref} 2> {output.log}; cat $(echo {output.log} | sed 's/.log//g')_*.fastq | gzip -3 > " + \
		"{output.rds}; cat $(echo {output.log} | sed 's/.log//g')_*.maf | gzip -3 > {output.maf}; rm $(echo {output.log} | sed " + \
		"'s/.log//g')_*.{{maf,ref,fastq}}"

rule simReads:
	input:
		ref = "../simulations/genomes/{genome}.fasta",
		model = "../software/pbsim2/data/{chem}.model"
	params:
		"{dRat}"
	output:
		rds = temp("../simulations/reads/{genome}_c{chem}_ep{dRat}.fastq"),
		log = "../simulations/reads/{genome}_c{chem}_ep{dRat}.log"
	shell:
		"pbsim --hmm_model {input.model} --difference-ratio {params} --prefix $(echo {output.log} | sed 's/.log//g') --depth " + \
		"10.0 {input.ref} 2> {output.log}; cat $(echo {output.log} | sed 's/.log//g')_*.fastq > {output.rds}; rm $(echo " + \
		"{output.log} | sed 's/.log//g')_*.{{ref,fastq}}"

rule simSeqs:
	output:
		"../simulations/simSeqs_n{n}_l{l}_m{m}_i{i}_d{d}.ssv"
	shell:
		"python3 scripts/SimSeqPairs.py -n {wildcards.n} -l {wildcards.l} -m {wildcards.m} -i {wildcards.i} -d {wildcards.d} > {output}"

rule measureSimilarity:
	input:
		"../simulations/simSeqs_n{simDesc}.ssv"
	params:
		simMes = "{mes}"
	output:
		"../simulations/scores_mes{mes}_n{simDesc}.txt"
	run:
		for p in open(input[0], 'r'):
			if not p.startswith("Random"):
				seqs = p.strip().split(' ')
				if params.simMes == "intersec":
					shell("src/CalcSim -i -a {seqs[0]} -b {seqs[1]} >> {output}")
				else:
					shell("src/CalcSim -l -a {seqs[0]} -b {seqs[1]} >> {output}")

rule simGeneSeqs:
	output:
		"../simulations/geneSeqs/simGeneSeqs_n{n}_l{l}_o{o}_m{m}_i{i}_d{d}.ssv"
	shell:
		"python3 scripts/SimSeqPairs.py -n {wildcards.n} -l {wildcards.l} -o {wildcards.o} -m {wildcards.m} -i {wildcards.i} -d " + \
		"{wildcards.d} > {output}"

rule genRandSeqs:
	output:
		"../simulations/randSeqs_n{n}_l{l}_m{m}_i{i}_d{d}.ssv"
	shell:
		"python3 scripts/SimSeqPairs.py -n {wildcards.n} -l {wildcards.l} -m {wildcards.m} -i {wildcards.i} -d {wildcards.d} > {output}"

rule constructSearchPairs:
	input:
		"../simulations/geneSeqs/simGeneSeqs_n{gn}_l{gl}_o{o}_m{m}_i{i}_d{d}.ssv",
		"../simulations/randSeqs_n{rn}_l{rl}_m{m}_i{i}_d{d}.ssv"
	output:
		"../simulations/searchPairs/searchPairs_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{i}_d{d}.txt"
	run:
		randSeqFile = open(input[1], 'r')
		randSeqFile.readline()

		for l in open(input[0], 'r'):
			if not l.startswith("Random"):
				genome = randSeqFile.readline().strip().split(' ')[0]
				isPattern = True

				for g in l.strip().split(' '):
					if isPattern:
						pattern = g
						isPattern = False
					else:
						genome += g + randSeqFile.readline().strip().split(' ')[0]

				shell("echo {pattern} {genome} >> {output}")

rule searchHomologies:
	input:
		"../simulations/searchPairs/searchPairs_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{i}_d{d}.txt"
	params:
		c = "{c}",
		u = "{u}",
		thres = "{t}"
	output:
		"../simulations/homologies/homologies_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{i}_d{d}_c{c}_u{u}_t{t}_pPy.txt"
	wildcard_constraints:
		t = "[0-9]+"
	run:
		i = 0

		for l in open(input[0], 'r'):
			i += 1
			shell("echo Pair {i} >> {output}")
			l = l.strip()
			pattern, text = l.split(' ')
			shell("python3 scripts/FindThoms.py -p {pattern} -s {text} -c {params.c} -u {params.u} -t {params.thres} >> {output}")

rule searchHomologiesIndvdlyWthScrptImpl:
	input:
		pttn = "../simulations/searchPairs/searchPairs_gn{desc}_qry{n}.fasta",
		txt = "../simulations/searchPairs/searchPairs_gn{desc}_ref{n}.fasta"
	params:
		c = "{c}",
		u = "{u}",
		thres = "{t}"
	output:
		"../simulations/homologies/homologies_gn{desc}_c{c}_u{u}_t{t}_n{n}_pPy.txt"
	shell:
		"python3 scripts/FindThoms.py -p {input.pttn} -s {input.txt} -c {params.c} -u {params.u} -t {params.thres} > {output}"

rule searchHomologiesWthCppImpl:
	input:
		pttn = "../simulations/searchPairs/searchPairs_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{i}_d{d}_qry{n}.fasta",
		txt = "../simulations/searchPairs/searchPairs_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{i}_d{d}_ref{n}.fasta"
	params:
		c = "{c}",
		u = "{u}",
		thres = "{t}"
	output:
		"../simulations/homologies/homologies_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{i}_d{d}_c{c}_u{u}_t{t}_n{n}_pCpp.txt"
	shell:
		"src/FindThoms -p {input.pttn} -s {input.txt} -c {params.c} -u {params.u} -t {params.thres} > {output}"

# rule aggregateResults:
# 	input:
# 		genAggList
# 	output:
# 		"../simulations/homologies/homologies_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{i}_d{d}_c{c}_u{u}_t{t}_pCpp.txt"
# 	run:
# 		files2Agg = sorted(glob('../simulations/searchPairs/searchPairs_gn100_rn400_gl1000_rl500_o3_m0_i0_d0_ref*.fasta'), key=\
# 			lambda n: int(n.split("ref")[1].split('.')[0]))
# 		shell("cat {' '.join(files2Agg)} > {output}")

rule createFASTApairs:
	input:
		"../simulations/searchPairs/searchPairs_{desc}.txt"
	params:
		"{i}"
	output:
		ref = "../simulations/searchPairs/searchPairs_{desc}_ref{i}.fasta",
		qry = "../simulations/searchPairs/searchPairs_{desc}_qry{i}.fasta"
	run:
		c = 0
		for l in open(input[0], 'r'):
			if c == int(params[0]):
				pattern, text = l.strip().split(' ')
				SeqIO.write([SeqRecord(Seq(text), id=f"SearchPair_{params[0]}_ref")], output[0], "fasta")
				SeqIO.write([SeqRecord(Seq(pattern), id=f"SearchPair_{params[0]}_qry")], output[1], "fasta")

			c += 1

rule runNucmer:
	input:
		ref = "../simulations/searchPairs/searchPairs_{desc}_d{d}_ref{i}.fasta",
		qry = "../simulations/searchPairs/searchPairs_{desc}_d{d}_qry{i}.fasta"
	output:
		"../simulations/nucmer/nucmerAlignments_{desc}_d{d, [0,1].[0-9]+}_p{i}.delta"
	shell:
		"/homes/tischulz/usr/local/bin/nucmer --prefix ../simulations/nucmer/nucmerAlignments_{wildcards.desc}_p{wildcards.i} " + \
		"{input.ref} {input.qry}"

rule runMaxmatchNucmer:
	input:
		ref = "../simulations/searchPairs/searchPairs_{desc}_ref{i}.fasta",
		qry = "../simulations/searchPairs/searchPairs_{desc}_qry{i}.fasta"
	output:
		"../simulations/nucmer/nucmerAlignments_{desc}_maxmatch_p{i}.delta"
	shell:
		"/homes/tischulz/usr/local/bin/nucmer --prefix ../simulations/nucmer/nucmerAlignments_{wildcards.desc}_maxmatch_p" + \
		"{wildcards.i} --maxmatch {input.ref} {input.qry}"

rule getAlignmentCoords:
	input:
		"../simulations/nucmer/nucmerAlignments_{desc}.delta"
	output:
		"../simulations/nucmer/nucmerAlignments_{desc}.coords"
	shell:
		"/homes/tischulz/usr/local/bin/show-coords {input} > {output}"
