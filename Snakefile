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

def genTextNames(wcs):
	divRates = config['divergenceRates']
	subIndelRats = config['substitutionIndelRates']
	textNames = []
	rl = 2 * config['templLen'] * (config['copyNumber'] + 1)

	for d in divRates:
		for s in subIndelRats:
			dsCombi = f"m{(d * s):.3f}_d{d}_i{d}"
			textNames += expand("../simulations/text_rl{l}_rid{i}_" + dsCombi + "_cn{C}.fasta", l=rl, i=range(100), C=\
				config['copyNumber'])

	return textNames

def genHomFiles(wcs):
	divRates = config['divergenceRates']
	subIndelRats = config['substitutionIndelRates']
	homFileNames = []
	rl = 2 * config['templLen'] * (config['copyNumber'] + 1)
	nbRds = int((config['readCoverage'] * config['templLen']) / config['readLength'])

	for d in divRates:
		for s in subIndelRats:
			dsCombi = f"m{(d * s):.3f}_d{d}_i{d}"
			homFileNames += expand("../simulations/homologies/homologies_rl{rl}_rid{i}_" + dsCombi + "_cn{C}_tl{tl}_chP6C4_" + \
				"ep0:0:0_ri{ri}_c1_u1_t0.txt", rl=rl, i=range(100), C=config['copyNumber'], tl=config['templLen'], ri=range(nbRds))

	return homFileNames

def genMinimap2Files(wcs):
	divRates = config['divergenceRates']
	subIndelRats = config['substitutionIndelRates']
	miniFileNames = []
	rl = 2 * config['templLen'] * (config['copyNumber'] + 1)
	nbRds = int((config['readCoverage'] * config['templLen']) / config['readLength'])

	for d in divRates:
		for s in subIndelRats:
			dsCombi = f"m{(d * s):.3f}_d{d}_i{d}"
			miniFileNames += expand("../simulations/minimap2Res/textmappings_rl{rl}_rid{i}_" + dsCombi + "_cn{C}_tl{tl}_chP6C4_" + \
				"ep0:0:0_ri{ri}.sam", rl=rl, i=range(100), C=config['copyNumber'], tl=config['templLen'], ri=range(nbRds))

	return miniFileNames

def genWinnowmap2Files(wcs):
	divRates = config['divergenceRates']
	subIndelRats = config['substitutionIndelRates']
	winFileNames = []
	rl = 2 * config['templLen'] * (config['copyNumber'] + 1)
	nbRds = int((config['readCoverage'] * config['templLen']) / config['readLength'])

	for d in divRates:
		for s in subIndelRats:
			dsCombi = f"m{(d * s):.3f}_d{d}_i{d}"
			winFileNames += expand("../simulations/Winnowmap2Res/textmappings_rl{rl}_rid{i}_" + dsCombi + "_cn{C}_tl{tl}_chP6C4_" + \
				"ep0:0:0_ri{ri}.sam", rl=rl, i=range(100), C=config['copyNumber'], tl=config['templLen'], ri=range(nbRds))

	return winFileNames

rule all:
	input:
		#Tests for DP script
		# expand("../simulations/homologies/homologies_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{m}_d{m}_{sp}_t0_pPy.txt", gn=\
		# 	config['nbSimSeqs'], rn=NB_RANDSEQS, gl=config['geneLen'], rl=config['randSeqLen'], o=config['nbCpys'], m=\
		# 	config['mutationRates'], sp=config['scoringPatterns']),
		# expand("../simulations/homologies/homologies_gn{gn}_rn200_gl{gl}_rl1000_o1_m{m}_i{m}_d{m}_{sp}_t0_pPy.txt", gn=\
		# 	config['nbSimSeqs'], gl=config['geneLen'], m=config['mutationRates'], sp=config['scoringPatterns']),
		#Nucmer runs
		# expand("../simulations/nucmer/nucmerAlignments_gn{gn}_rn200_gl{gl}_rl1000_o1_m{m}_i{m}_d{m}_p{i}.coords", \
		# 	gn=config['nbSimSeqs'], gl=config['geneLen'], m=config['mutationRates'], i=range(config['nbSimSeqs'])),
		# expand("../simulations/nucmer/nucmerAlignments_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{m}_d{m}_p{i}.coords", gn=\
		# 	config['nbSimSeqs'], rn=NB_RANDSEQS, gl=config['geneLen'], rl=config['randSeqLen'], o=config['nbCpys'], m=\
		# 	config['mutationRates'], i=range(config['nbSimSeqs'])),
		# expand("../simulations/nucmer/nucmerAlignments_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{m}_d{m}_maxmatch_p{i}.coords", gn=\
		# 	config['nbSimSeqs'], rn=NB_RANDSEQS, gl=config['geneLen'], rl=config['randSeqLen'], o=config['nbCpys'], m=\
		# 	config['mutationRates'], i=range(config['nbSimSeqs'])),
		#Tests for DP C++ implementation
		# expand("../simulations/homologies/homologies_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{m}_d{m}_{sp}_t0_n{n}_p{p}.uni", gn=\
		# 	config['nbSimSeqs'], rn=NB_RANDSEQS, gl=config['geneLen'], rl=config['randSeqLen'], o=config['nbCpys'], m=\
		# 	config['mutationRates'], sp=config['scoringPatterns'], n=range(100), p=config['implType']),
		# expand("../simulations/homologies/homologies_gn{gn}_rn200_gl{gl}_rl1000_o1_m{m}_i{m}_d{m}_{sp}_t0_n{n}_p{p}.uni", gn=\
		# 	config['nbSimSeqs'], gl=config['geneLen'], m=config['mutationRates'], sp=config['scoringPatterns'], n=range(100), p=\
		# 	config['implType']),
		#minimap2 and bwamem runs
		# matchProgCallToSeed,
		#
		genHomFiles,
		genMinimap2Files,
		genWinnowmap2Files

rule convertMultiFastq2SinglFastq:
	input:
		"{something}.fastq"
	output:
		expand("{s}_{i}.fastq", s="{something}", i=range(int((config['readCoverage'] * config['templLen']) / config['readLength'])))
	shell:
		"python3 scripts/convMulFq2SinFq.py {input}"

rule convertMultiFastq2SinglFasta:
	input:
		"{something}.fastq"
	output:
		expand("{s}_{i}.fasta", s="{something}", i=range(int((config['readCoverage'] * config['templLen']) / config['readLength'])))
	shell:
		"python3 scripts/convMulFq2SinFa.py {input}"

rule assembleText:
	input:
		randSeq = "../simulations/randSeq_l{rl}_rid{i}.fasta",
		cpys = expand("../simulations/randSeqCopy_l{tl}_rid{i}_m{m}_d{d}_i{l}_cn{c}.fasta", tl=config['templLen'], i="{i}", m="{m}"\
			, d="{d}", l="{l}", c=range(int(config['copyNumber'])))
	output:
		"../simulations/text_rl{rl}_rid{i}_m{m}_d{d}_i{l}_cn{c}.fasta"
	shell:
		"python3 scripts/assembleText.py -b {input.randSeq} -i {input.cpys} -o {output}"

rule mutateTemplate:
	input:
		"../simulations/randSeq_{desc}.fasta"
	params:
		subR = "{m}",
		delR = "{d}",
		iLen = "{l}",
		rid = "{i}"
	output:
		"../simulations/randSeqCopy_{desc}_m{m}_d{d}_i{l}_cn{i}.fasta"
	shell:
		"python3 scripts/MutateSeq.py -m {params.subR} -d {params.delR} -i {params.iLen} -t {input} -o {output}"

rule genRandTemplate:
	params:
		length = "{l}",
		repId = "{i}"
	output:
		"../simulations/randSeq_l{l}_rid{i}.fasta"
	shell:
		"python3 scripts/GenRandSeq.py -l {params.length} -o {output}"

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

rule runWinnowmap2:
	input:
		ref = "../simulations/text_rl{rl}_rid{ri}_m{desc}.fasta",
		qry = "../simulations/reads/reads_tl{tl}_rid{ri}_ch{rdDesc}_{i}.fastq",
		cnts = "../simulations/repKmers_k15_rl{rl}_rid{ri}_m{desc}.txt"
	output:
		res = "../simulations/Winnowmap2Res/textmappings_rl{rl}_rid{ri}_m{desc}_tl{tl}_ch{rdDesc}_ri{i}.sam",
		bench = "../benchmarks/benchWinnowmap2_rl{rl}_rid{ri}_m{desc}_tl{tl}_ch{rdDesc}_ri{i}.txt"
	shell:
		"/usr/bin/time -v ../software/Winnowmap/bin/winnowmap -W {input.cnts} -ax map-pb {input.ref} {input.qry} > {output.res}" + \
		" 2> {output.bench}"

rule printCounts:
	input:
		"../simulations/merylDB_k{k}_rl{desc}"
	output:
		temp("../simulations/repKmers_k{k}_rl{desc}.txt")
	shell:
		"../software/Winnowmap/bin/meryl print greater-than distinct=0.9998 {input} > {output}"

rule countKmers:
	input:
		"../simulations/text_rl{desc}.fasta"
	params:
		"{k}"
	output:
		temp(directory("../simulations/merylDB_k{k}_rl{desc}"))
	shell:
		"../software/Winnowmap/bin/meryl count k={params} output {output} {input}"

rule runMinimap2onPacBio:
	input:
		ref = "../simulations/text_rl{rl}_rid{ri}_m{desc}.fasta",
		qry = "../simulations/reads/reads_tl{tl}_rid{ri}_ch{rdDesc}_{i}.fastq"
	output:
		res = "../simulations/minimap2Res/textmappings_rl{rl}_rid{ri}_m{desc}_tl{tl}_ch{rdDesc}_ri{i}.sam",
		bench = "../benchmarks/benchMinimap2_rl{rl}_rid{ri}_m{desc}_tl{tl}_ch{rdDesc}_ri{i}.txt"
	shell:
		"/usr/bin/time -v minimap2 -ax map-hifi {input.ref} {input.qry} > {output.res} 2> {output.bench}"

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

rule simPattern:
	input:
		ref = "../simulations/randSeq_l{desc}.fasta",
		model = "../software/pbsim2/data/{chem}.model"
	params:
		erR = "{erR}",
		dep = config['readCoverage'],
		rdLen = config['readLength']
	output:
		rds = "../simulations/reads/reads_tl{desc}_ch{chem}_ep{erR}.fastq",
		log = "../simulations/reads/reads_tl{desc}_ch{chem}_ep{erR}.log"
	wildcard_constraints:
		erR="[0-9]+:[0-9]+:[0-9]+"
	shell:
		"pbsim --hmm_model {input.model} --difference-ratio {params} --prefix $(echo {output.log} | sed 's/.log//g') --depth " + \
		"{params.dep} --length-min {params.rdLen} --length-max {params.rdLen} {input.ref} 2> {output.log}; cat $(echo " + \
		"{output.log} | sed 's/.log//g')_*.fastq > {output.rds}; rm $(echo {output.log} | sed 's/.log//g')_*.{{ref,fastq}}"

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
		pttn = "../simulations/reads/reads_tl{tl}_rid{ri}_ch{chem}_ep{erR}_{i}.fasta",
		txt = "../simulations/text_rl{rl}_rid{ri}_m{desc}.fasta"
	params:
		c = "{c}",
		u = "{u}",
		thres = "{t}",
		k = config['k']
	output:
		homs = "../simulations/homologies/homologies_rl{rl}_rid{ri}_m{desc}_tl{tl}_ch{chem}_ep{erR}_ri{i}_c{c}_u{u}_t{t}.txt",
		bench = "../benchmarks/benchFindThoms_rl{rl}_rid{ri}_m{desc}_tl{tl}_ch{chem}_ep{erR}_ri{i}_c{c}_u{u}_t{t}.txt"
	shell:
		"/usr/bin/time -v src/FindThoms -p {input.pttn} -s {input.txt} -k {params.k} -c {params.c} -u {params.u} -t {params.thres} > " + \
		"{output.homs} 2> {output.bench}"

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
