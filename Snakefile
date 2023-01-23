configfile: 'config.yaml'

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from os.path import exists
from random import randrange, seed
from sys import maxsize
from glob import glob

NB_RANDSEQS = (config['nbCpys'] + 1) * config['nbSimSeqs']
# T2T_READ_LEN = [len(r.seq) for r in SeqIO.parse(open(f"../simulations/reads/t2thumanChrY_cP6C4_ep6:50:54_s322235950486831966.fasta", 'r'), "fasta")]
DIV_RATE = 0.01
SUB_RATE = DIV_RATE / 3
INS_RATE = DIV_RATE / 3
DEL_RATE = DIV_RATE / 3
SUB_ERR = 0.002 * 6 / 110.
INS_ERR = 0.002 * 50 / 110.
DEL_ERR = 0.002 * 54 / 110.

#Initialize random number generator with global seed for this workflow
seed(config['globalSeed'])

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
				"ep0:0:0_ri{ri}_k{k}_c1_u1_t0.txt", rl=rl, i=range(100), C=config['copyNumber'], tl=config['templLen'], ri=\
				range(nbRds), k=config['k'])
			homFileNames += expand("../simulations/homologies/homologies_rl{rl}_rid{i}_" + dsCombi + "_cn{C}_tl{tl}_chP6C4_" + \
				"ep0:0:0_ri{ri}_k15_c1_u1_t-19999.txt", rl=rl, i=range(100), C=config['copyNumber'], tl=config['templLen'], ri=\
				range(nbRds))

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
				"ep0:0:0_k{k}_ri{ri}.sam", rl=rl, i=range(100), C=config['copyNumber'], tl=config['templLen'], k=config['k'], ri=\
				range(nbRds))

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
				"ep0:0:0_k{k}_ri{ri}.sam", rl=rl, i=range(100), C=config['copyNumber'], tl=config['templLen'], k=config['k'], ri=\
				range(nbRds))

	return winFileNames

def genSampleFileNames(wcs):
	scoreFiles = []
	sr = SUB_RATE
	ir = INS_RATE
	dr = DEL_RATE
	ser = SUB_ERR
	ier = INS_ERR
	der = DEL_ERR

	for l in config['readLens']:
		for i in range(config['randomSampleSize']):
			seed = randrange(maxsize)
			scoreFiles.append(f"../simulations/expValExp/scores/readRefGlobIntSecScore_se{seed}_sl{l}_sr{sr}_ir{ir}_dr{dr}" + \
				f"_ser{ser}_ier{ier}_der{der}_k15_r0.1_c1_u1.txt")

	tl = config['templLen']
	k = config['minimap2DefaultK']

	# for i in range(config['randomSampleSize']):
	# 	for d in config['divergenceRates']:
	# 		for s in config['substitutionIndelRates']:
	# 			sd = randrange(maxsize)
	# 			sr = f"{(d * s):.3f}"
	# 			scoreFiles.append(f"../simulations/expValExp/scores/noDupGlobIntSecScore_se{sd}_sl{tl}_sr{sr}_ir{d}_dr{d}_k{k}" + \
	# 				"_r0.1_c1_u1.txt")

	# for i in range(config['randomSampleSize']):
	# 	for d in config['divergenceRates']:
	# 		for s in config['substitutionIndelRates']:
	# 			sd = randrange(maxsize)
	# 			sr = f"{(d * s):.3f}"
	# 			scoreFiles.append(f"../simulations/expValExp/scores/dupGlobIntSecScore_se{sd}_sl{tl}_sr{sr}_ir{d}_dr{d}_k{k}" + \
	# 				"_r0.1_c1_u1.txt")

	# for i in range(config['randomSampleSize']):
	# 	for d in config['divergenceRates']:
	# 		for s in config['substitutionIndelRates']:
	# 			sd = randrange(maxsize)
	# 			sr = f"{(d * s):.3f}"
	# 			scoreFiles.append(f"../simulations/expValExp/scores/noDupGlobIntSecScore_se{sd}_sl{config['readLength']}_sr{sr}" + \
	# 				f"_ir{d}_dr{d}_k{k}_r0.1_c1_u1.txt")
	# 			scoreFiles.append(f"../simulations/expValExp/scores/dupGlobIntSecScore_se{sd}_sl{config['readLength']}_sr{sr}" + \
	# 				f"_ir{d}_dr{d}_k{k}_r0.1_c1_u1.txt")

	return scoreFiles

def genReadScoreFiles(wcs):
	scrFiles = []

	for l in config['readLens']:
		for i in range(config['randomSampleSize']):
			sd = randrange(maxsize)
			scrFiles.append(f"../simulations/reads/randSeq_l{int(l * config['refRdRatio'])}_rid{i}" + \
				f"_cP6C4_ep6:50:54_s{sd}_d1_l{l}-{l}.maf.gz")
			# scrFiles.append(f"../simulations/expValExp/scores/randSeq_l{int(l * config['refRdRatio'])}_rid{i}GlobIntSecScore_se" + \
			# 	f"{sd}_cP6C4_ep6:50:54_d1_l{l}-{l}_c1_u1.txt")

	return scrFiles

def genReadScoreFiles2(wcs):
	scrFiles = []

	for l in [1600, 3200, 6400]:
		for i in range(config['randomSampleSize']):
			sd = randrange(maxsize)
			scrFiles.append(f"../simulations/expValExp/scores/randSeq_l{l}_rid{i}GlobIntSecScore_se{sd}_cP6C4_ep6:50:54_d1_l{l}-" + \
				f"{l}_c1_u1.txt")
	return scrFiles

def genParaSailResFiles(wcs):
	rdLen = T2T_READ_LEN[int(wcs.rdId) - 1]
	refLen = config['t2tChrYlen']
	maxPieceSize = 2000000
	res = []

	for ps in [s for s in range(0, refLen, 2000000 - 2 * rdLen) if s + 2 * rdLen < refLen]:
		pe = ps + min(ps + maxPieceSize, refLen)

		res.append(f"../simulations/parasailMappings/{wcs.gn}_ra{ps}-{pe}_c{wcs.ch}_ep{wcs.eP}_s{wcs.sd}_ri{wcs.rdId}.pr")

	return res	

READ_SEED = randrange(maxsize)

rule all:
	input:
		#Expectation value estimation
		# genSampleFileNames,
		# genReadScoreFiles,
		# expand("../simulations/expValExp/dists/globEditDistance_srandSeq_l{l}_rid{n}_trandSeq_l{l}_rid{n}_m{m}_d{d}_i{i}_cn0_m" + \
		# 	"{se}_d{de}_i{ie}_cn0.txt", l=config['readLens'], n=range(config['randomSampleSize']), m=SUB_RATE, d=DEL_RATE, i=\
		# 	INS_RATE, se=SUB_ERR, de=DEL_ERR, ie=INS_ERR),
		#Testing
		# genReadScoreFiles2,
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
		#DP implementation benchmark on human data
		# expand("../benchmarks/benchFindThoms_humanChr20_refined_onlyCapitalNucs_ep{e}_s1657921994_rr{i}_k15_c1_u1_t-1000.txt", e=\
		# config['errorPatterns'], i=range(10)),
		#Benchmark on real data
		# f"../simulations/reads/t2thumanChrY_sr{SUB_ERR}_dr{DEL_ERR}_i{INS_ERR}_sd{randrange(maxsize)}_lmn" + \
		# f"{config['pbsimLenMin']}_lmx{config['pbsimLenMax']}_lavg{config['pbsimLenAvg']}_ls{config['pbsimLenStd']}_dp10.fasta",
		# expand("../simulations/edlibMappings/t2thumanChrY_sr{sr}_dr{dr}_i{ir}_sd{s}_lmn" + \
		# "{mn}_lmx{mx}_lavg{m}_ls{sd}_dp10_ri{i}.er", \
		# sr=SUB_ERR, dr=DEL_ERR, ir=INS_ERR, s=READ_SEED, mn=config['pbsimLenMin'], mx=config['pbsimLenMax'], m=\
		# config['pbsimLenAvg'], sd=config['pbsimLenStd'], i=range(69401)),
		f"../simulations/homologies/homologies_t2thumanChrY_sr{SUB_ERR}_dr{DEL_ERR}_i{INS_ERR}_sd{READ_SEED}_lmn" + \
		f"{config['pbsimLenMin']}_lmx{config['pbsimLenMax']}_lavg{config['pbsimLenAvg']}_ls{config['pbsimLenStd']}_dp10_k15_" + \
			f"hr{config['hashRate']}_c1_u1_de{0.05226723}_in{-116.02672267226808}.txt",
		# expand("../benchmarks/benchFindThoms_t2thumanChrY_sr{sr}_dr{dr}_i{ie}_sd{sd}_lmn{mn}_lmx{mx}_lavg{m}_ls{s}_dp10_k15_" + \
		# 	"hr{hr}_c1_u1_d{de}_in{it}_rep{i}.txt", sr=SUB_ERR, dr=DEL_ERR, ie=INS_ERR, sd=READ_SEED, mn=config['pbsimLenMin'], mx=\
		# 	config['pbsimLenMax'], m=config['pbsimLenAvg'], s=config['pbsimLenStd'], hr=config['hashRate'], de=0.05226723, it=\
		# 	-116.02672267226808, i=range(config['benchRepRuns'])),
		# f"../simulations/minimap2Res/t2thumanChrY_sr{SUB_ERR}_dr{DEL_ERR}_i{INS_ERR}_sd{READ_SEED}_lmn{config['pbsimLenMin']}" + \
		# f"_lmx{config['pbsimLenMax']}_lavg{config['pbsimLenAvg']}_ls{config['pbsimLenStd']}_dp10_k15.sam.gz",
		# expand("../benchmarks/benchMinimap2_t2thumanChrY_sr{sr}_dr{dr}_i{ie}_sd{sd}_lmn{mn}_lmx{mx}_lavg{m}_ls{s}_dp10_k15_rep" + \
		# 	"{i}.txt", sr=SUB_ERR, dr=DEL_ERR, ie=INS_ERR, sd=READ_SEED, mn=config['pbsimLenMin'], mx=config['pbsimLenMax'], m=\
		# 	config['pbsimLenAvg'], s=config['pbsimLenStd'], i=range(config['benchRepRuns'])),
		# f"../simulations/Winnowmap2Res/t2thumanChrY_sr{SUB_ERR}_dr{DEL_ERR}_i{INS_ERR}_sd{READ_SEED}_lmn{config['pbsimLenMin']}" + \
		# f"_lmx{config['pbsimLenMax']}_lavg{config['pbsimLenAvg']}_ls{config['pbsimLenStd']}_dp10_k15.sam.gz",
		# expand("../benchmarks/benchWinnowmap2_t2thumanChrY_sr{sr}_dr{dr}_i{ie}_sd{sd}_lmn{mn}_lmx{mx}_lavg{m}_ls{s}_dp10_k15" + \
		# 	"_rep{i}.txt", sr=SUB_ERR, dr=DEL_ERR, ie=INS_ERR, sd=READ_SEED, mn=config['pbsimLenMin'], mx=config['pbsimLenMax'], m=\
		# 	config['pbsimLenAvg'], s=config['pbsimLenStd'], i=range(config['benchRepRuns'])),
		# expand("../benchmarks/benchFindThoms_t2thumanChrY_chP6C4_ep6:50:54_s322235950486831966_k15_hr0.2_c1_u1_t0_rep{i}.txt", i=\
		# 	range(config['benchRepRuns'])),
		# expand("../benchmarks/benchMinimap2_t2thumanChrY_cP6C4_ep6:50:54_s322235950486831966_k15_rep{i}.txt", i=\
		# 	range(config['benchRepRuns'])),
		# expand("../benchmarks/benchWinnowmap2_t2thumanChrY_cP6C4_ep6:50:54_s322235950486831966_k15_rep{i}.txt", i=\
		# 	range(config['benchRepRuns'])),
		#
		# genHomFiles,
		# genMinimap2Files,
		# genWinnowmap2Files

rule calcGlobEditDist:
	input:
		s = "../simulations/{sNm}.fasta",
		t = "../simulations/{tNm}.fasta"
	output:
		"../simulations/expValExp/dists/globEditDistance_s{sNm}_t{tNm}.txt"
	shell:
		"FindSimSeqs/CalcGlobEditDistance {input.s} {input.t} > {output}"

rule aggregateParasailRes:
	input:
		genParaSailResFiles
	output:
		"../simulations/parasailMappings/{gn}_c{ch}_ep{eP}_s{sd}_ri{rdId}.tpr"
	shell:
	 	"python3 scripts/aggrParasailRes.py {input} > {output}"

rule runParasail:
	input:
		ref = "../simulations/genomes/{genome}_ra{range}.fasta",
		qry = "../simulations/reads/{genome}_c{rdDesc}_ri{rdId}.fasta"
	output:
		"../simulations/parasailMappings/{genome}_ra{range}_c{rdDesc}_ri{rdId}.pr"
	shell:
		"FindSimSeqs/FindSimSeqs {input.qry} {input.ref} > {output}"

rule splitRef:
	input:
		"../simulations/genomes/{genome}.fasta"
	params:
		start = "{s}",
		end = "{e}"
	output:
		"../simulations/genomes/{genome}_ra{s}-{e}.fasta"
	shell:
		"python3 scripts/getSubstring.py -i {input} -s {params.start} -e {params.end} -o {output}"

rule runEdlib:
	input:
		ref = "../simulations/genomes/{genome}.fasta",
		qry = "../simulations/reads/{genome}_sr{rdDesc}_ri{rdId}.fasta"
	output:
		"../simulations/edlibMappings/{genome}_sr{rdDesc}_ri{rdId}.er"
	shell:
		"FindSimSeqs/FindSimSeqs {input.qry} {input.ref} > {output}"

rule divideReads:
	input:
		"../simulations/reads/{rdFileName}.fasta"
	params:
		"{rdId}"
	output:
		temp("../simulations/reads/{rdFileName}_ri{rdId}.fasta")
	shell:
		"python3 scripts/getSeq.py -s {input} -i {params} -o {output}"

rule calculateGlobalIntersectionSimilarity:
	input:
		"../simulations/expValExp/sketches/{dupInfo}KskSeqPairs_se{desc}.sk"
	output:
		"../simulations/expValExp/scores/{dupInfo}GlobIntSecScore_se{desc}_c1_u1.txt"
	shell:
		"python3 scripts/CalcGlobInterSim.py -p {input} > {output}"

rule generateReadRefSketchPairs:
	params:
		sd = "{sd}",
		sl = "{sl}",
		sr = "{sr}",
		ir = "{ir}",
		dr = "{dr}",
		ser = "{ser}",
		ier = "{ier}",
		der = "{der}",
		k = "{k}",
		r = "{hRat}"
	output:
		"../simulations/expValExp/sketches/readRefKskSeqPairs_se{sd}_sl{sl}_sr{sr}_ir{ir}_dr{dr}_ser{ser}" + \
		"_ier{ier}_der{der}_k{k}_r{hRat}.sk"
	shell:
		"python3 scripts/genRdRfSks.py -s {params.sd} -l {params.sl} -m {params.sr} -i {params.ir} -d {params.dr} -se " + \
		"{params.ser} -ie {params.ier} -de {params.der} -k {params.k} -r {params.r} > {output}"

rule generateKmerSketchPairs:
	output:
		"../simulations/expValExp/sketches/dupKskSeqPairs_se{seed}_sl{seqLen}_sr{subRate}_ir{insRate}_dr{delRate}_k{k}_r{hRat}.sk"
	shell:
		"python3 scripts/GenKsketchPair.py -s {wildcards.seed} -l {wildcards.seqLen} -m {wildcards.subRate} -i " + \
		"{wildcards.insRate} -d {wildcards.delRate} -k {wildcards.k} -r {wildcards.hRat} > {output}"

rule generateNoDuplicateKmerSketchPairs:
	output:
		"../simulations/expValExp/sketches/noDupKskSeqPairs_se{seed}_sl{seqLen}_sr{subRate}_ir{insRate}_dr{delRate}_k{k}_r{hRat}.sk"
	shell:
		"python3 scripts/GenKsketchPair.py -n -s {wildcards.seed} -l {wildcards.seqLen} -m {wildcards.subRate} -i " + \
		"{wildcards.insRate} -d {wildcards.delRate} -k {wildcards.k} -r {wildcards.hRat} > {output}"

rule generateReadSketchPairs:
	input:
		tmpl = "../simulations/genomes/{genome}.fasta",
		maf = "../simulations/reads/{genome}_c{desc}_s{sd}_d{desc1}.maf.gz"
	output:
		"../simulations/expValExp/sketches/{genome}KskSeqPairs_se{sd}_c{desc}_d{desc1}.sk"
	shell:
		"python3 scripts/genRdSkPrs.py {input.tmpl} {input.maf} > {output}"

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
		randSeq = "../simulations/genomes/randSeq_l{rl}_rid{i}.fasta",
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
		"../simulations/randSeq_{desc}_m{m}_d{d}_i{l}_cn{i}.fasta"
	shell:
		"python3 scripts/MutateSeq.py -m {params.subR} -d {params.delR} -i {params.iLen} -t {input} -o {output}"

rule genRandTemplate:
	params:
		length = "{l}",
		repId = "{i}"
	output:
		"../simulations/randSeq_l{l}_rid{i}.fasta"
	wildcard_constraints:
		i="[0-9]+"
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

rule saveWinnowmap2Result:
	input:
		"../simulations/Winnowmap2Res/{genome}_sr{desc}_k{k}_rep0.sam.gz"
	output:
		"../simulations/Winnowmap2Res/{genome}_sr{desc}_k{k}.sam.gz"
	wildcard_constraints:
		k = "[0-9]+"
	shell:
		"mv {input} {output}"

rule runWinnowmap2onRealGenomeFASTA:
	input:
		ref = "../simulations/genomes/{genome}.fasta",
		qry = "../simulations/reads/{genome}_sr{desc}.fasta",
		cnts = "../simulations/repKmers_k{k}_{genome}.txt"
	params:
		k = "{k}",
		r = "{r}"
	output:
		res = temp("../simulations/Winnowmap2Res/{genome}_sr{desc}_k{k}_rep{r}.sam.gz"),
		bench = "../benchmarks/benchWinnowmap2_{genome}_sr{desc}_k{k}_rep{r}.txt"
	wildcard_constraints:
		r = "[0-9]+"
	shell:
		"/usr/bin/time -v ../software/Winnowmap/bin/winnowmap -W {input.cnts} -ax map-pb -k {params.k} {input.ref} {input.qry} " + \
		"2> {output.bench} | gzip -3 > {output.res}"

rule runWinnowmap2onRealGenome:
	input:
		ref = "../simulations/genomes/{genome}.fasta",
		qry = "../simulations/reads/{genome}_c{chem}_ep{dRat}_s{seed}.fastq.gz",
		cnts = "../simulations/repKmers_k{k}_{genome}.txt"
	params:
		k = "{k}",
		r = "{r}"
	output:
		res = temp("../simulations/Winnowmap2Res/{genome}_c{chem}_ep{dRat}_s{seed}_k{k}_rep{r}.sam.gz"),
		bench = "../benchmarks/benchWinnowmap2_{genome}_c{chem}_ep{dRat}_s{seed}_k{k}_rep{r}.txt"
	shell:
		"/usr/bin/time -v ../software/Winnowmap/bin/winnowmap -W {input.cnts} -ax map-pb -k {params.k} {input.ref} {input.qry} " + \
		"2> {output.bench} | gzip -3 > {output.res}"

rule runWinnowmap2:
	input:
		ref = "../simulations/text_rl{rl}_rid{ri}_m{desc}.fasta",
		qry = "../simulations/reads/reads_tl{tl}_rid{ri}_ch{rdDesc}_{i}.fastq",
		cnts = "../simulations/repKmers_k{k}_rl{rl}_rid{ri}_m{desc}.txt"
	params:
		"{k}"
	output:
		res = "../simulations/Winnowmap2Res/textmappings_rl{rl}_rid{ri}_m{desc}_tl{tl}_ch{rdDesc}_k{k}_ri{i}.sam",
		bench = "../benchmarks/benchWinnowmap2_rl{rl}_rid{ri}_m{desc}_tl{tl}_ch{rdDesc}_{k}_ri{i}.txt"
	shell:
		"/usr/bin/time -v ../software/Winnowmap/bin/winnowmap -W {input.cnts} -ax map-pb -k {params} {input.ref} {input.qry} > " + \
		"{output.res} 2> {output.bench}"

rule printCounts:
	input:
		"../simulations/merylDB_k{k}_{desc}"
	output:
		temp("../simulations/repKmers_k{k}_{desc}.txt")
	shell:
		"../software/Winnowmap/bin/meryl print greater-than distinct=0.9998 {input} > {output}"

rule countGenomeKmers:
	input:
		"../simulations/genomes/{genome}.fasta"
	params:
		"{k}"
	output:
		temp(directory("../simulations/merylDB_k{k}_{genome}"))
	shell:
		"../software/Winnowmap/bin/meryl count k={params} output {output} {input}"

rule countKmers:
	input:
		"../simulations/text_rl{desc}.fasta"
	params:
		"{k}"
	output:
		temp(directory("../simulations/merylDB_k{k}_rl{desc}"))
	shell:
		"../software/Winnowmap/bin/meryl count k={params} output {output} {input}"

rule saveMinimap2Result:
	input:
		"../simulations/minimap2Res/{genome}_sr{desc}_k{k}_rep0.sam.gz"
	output:
		"../simulations/minimap2Res/{genome}_sr{desc}_k{k}.sam.gz"
	wildcard_constraints:
		k = "[0-9]+"
	shell:
		"mv {input} {output}"

rule runMinimap2onRealGenomePacBioFASTA:
	input:
		ref = "../simulations/genomes/{genome}.fasta",
		qry = "../simulations/reads/{genome}_sr{desc}.fasta"
	params:
		k = "{k}",
		r = "{r}"
	output:
		res = temp("../simulations/minimap2Res/{genome}_sr{desc}_k{k}_rep{r}.sam.gz"),
		bench = "../benchmarks/benchMinimap2_{genome}_sr{desc}_k{k}_rep{r}.txt"
	wildcard_constraints:
		r = "[0-9]+"
	shell:
		"/usr/bin/time -v minimap2 -ax map-hifi -k {params.k} {input.ref} {input.qry} 2> {output.bench} | gzip -3 > {output.res}"

rule runMinimap2onRealGenomePacBio:
	input:
		ref = "../simulations/genomes/{genome}.fasta",
		qry = "../simulations/reads/{genome}_c{chem}_ep{dRat}_s{seed}.fastq.gz"
	params:
		k = "{k}",
		r = "{r}"
	output:
		res = temp("../simulations/minimap2Res/{genome}_c{chem}_ep{dRat}_s{seed}_k{k}_rep{r}.sam.gz"),
		bench = "../benchmarks/benchMinimap2_{genome}_c{chem}_ep{dRat}_s{seed}_k{k}_rep{r}.txt"
	shell:
		"/usr/bin/time -v minimap2 -ax map-hifi -k {params.k} {input.ref} {input.qry} 2> {output.bench} | gzip -3 > {output.res}"

rule runMinimap2onPacBio:
	input:
		ref = "../simulations/text_rl{rl}_rid{ri}_m{desc}.fasta",
		qry = "../simulations/reads/reads_tl{tl}_rid{ri}_ch{rdDesc}_{i}.fastq"
	params:
		"{k}"
	output:
		res = "../simulations/minimap2Res/textmappings_rl{rl}_rid{ri}_m{desc}_tl{tl}_ch{rdDesc}_k{k}_ri{i}.sam",
		bench = "../benchmarks/benchMinimap2_rl{rl}_rid{ri}_m{desc}_tl{tl}_ch{rdDesc}_k{k}_ri{i}.txt"
	shell:
		"/usr/bin/time -v minimap2 -ax map-hifi -k {params} {input.ref} {input.qry} > {output.res} 2> {output.bench}"

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
	wildcard_constraints:
		seed="[0-9]+"
	shell:
		"pbsim --hmm_model {input.model} --difference-ratio {params.dr} --prefix $(echo {output.log} | sed 's/.log//g') --depth" + \
		" 10.0 --seed {params.sd} {input.ref} 2> {output.log}; cat $(echo {output.log} | sed 's/.log//g')_*.fastq | gzip -3 > " + \
		"{output.rds}; cat $(echo {output.log} | sed 's/.log//g')_*.maf | gzip -3 > {output.maf}; rm $(echo {output.log} | sed " + \
		"'s/.log//g')_*.{{maf,ref,fastq}}"

rule simPattern:
	input:
		ref = "../simulations/genomes/randSeq_l{desc}.fasta",
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

rule simFixedLengthReads:
	input:
		ref = "../simulations/genomes/{genome}.fasta",
		model = "../software/pbsim2/data/{chem}.model"
	params:
		dr = "{dRat}",
		dp = "{dpth}",
		lmin = "{minl}",
		lmax = "{maxl}",
		sd = "{seed}",
	output:
		rds = temp("../simulations/reads/{genome}_c{chem}_ep{dRat}_s{seed}_d{dpth}_l{minl}-{maxl}.fastq.gz"),
		maf = temp("../simulations/reads/{genome}_c{chem}_ep{dRat}_s{seed}_d{dpth}_l{minl}-{maxl}.maf.gz"),
		log = "../simulations/reads/{genome}_c{chem}_ep{dRat}_s{seed}_d{dpth}_l{minl}-{maxl}.log"
	shell:
		"pbsim --hmm_model {input.model} --difference-ratio {params.dr} --prefix $(echo {output.log} | sed 's/.log//g') " + \
		"--depth {params.dp} --seed {params.sd} --length-min {params.lmin} --length-mean {params.lmin} --length-sd 0 " + \
		"--length-max {params.lmax} {input.ref} 2> {output.log}; cat $(echo {output.log} | sed 's/.log//g')_*.fastq | gzip -3 >" + \
		" {output.rds}; cat $(echo {output.log} | sed 's/.log//g')_*.maf | gzip -3 > {output.maf}; rm $(echo {output.log} | sed" + \
		" 's/.log//g')_*.{{maf,ref,fastq}}"

rule simReadsOwnScript:
	input:
		ref = "../simulations/genomes/{genome}.fasta"
	params:
		dp = "{dp}",
		lMin = "{lMin}",
		lMax = "{lMax}",
		lMean = "{lMn}",
		lStd = "{lStd}",
		subR = "{subR}",
		delR = "{delR}",
		insR = "{insR}",
		sd = "{sd}"
	output:
		rds = "../simulations/reads/{genome}_sr{subR}_dr{delR}_i{insR}_sd{sd}_lmn{lMin}_lmx{lMax}_lavg{lMn}_ls{lStd}_dp{dp}.fasta"
	wildcard_constraints:
		dp = "[0-9]+"
	shell:
		"python3 scripts/simReads.py -dp {params.dp} -lmn {params.lMin} -lmx {params.lMax} -lavg {params.lMean} -ls " + \
		"{params.lStd} -r {input} -sr {params.subR} -dr {params.delR} -ir {params.insR} -sd {params.sd} -o {output}"

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

rule drawRandomRead:
	input:
		"../simulations/reads/humanChr20_cR103_ep{errorPattern}_s{seed}.fastq.gz"
	params:
		"{i}"
	output:
		"../simulations/reads/humanChr20_cR103_ep{errorPattern}_s{seed}_rr{i}.fasta"
	shell:
		"SCRIPT_INPUT=$(echo {input} | sed 's/.gz//g'); gunzip -c {input} > $SCRIPT_INPUT; python3 scripts/DrawRandSeqRec.py " + \
		"$SCRIPT_INPUT {output}"

rule convertCompressedFastq2Fasta:
	input:
		"{pth}.fastq.gz"
	output:
		"{pth}.fasta"
	run:
		fqName = input[0].replace(".gz", "")
		shell("gunzip -c {input} > %s" %fqName)
		SeqIO.write(SeqIO.parse(fqName, "fastq"), open(output[0], 'w'), "fasta")
		shell("rm %s" %fqName)

rule saveFindThomsResult:
	input:
		"../simulations/homologies/homologies_{genome}_{desc}_k{k}_hr{hr}_c{c}_u{u}_de{d}_in{i}_rep0.txt"
	output:
		"../simulations/homologies/homologies_{genome}_{desc}_k{k}_hr{hr}_c{c}_u{u}_de{d}_in{i}.txt"
	wildcard_constraints:
		i = "-?[0-9]+\.?[0-9]*"
	shell:
		"mv {input} {output}"

rule searchReadHomologies:
	input:
		rds = "../simulations/reads/{genome}_{desc}.fasta",
		txt = "../simulations/genomes/{genome}.fasta"
	params:
		c = "{c}",
		u = "{u}",
		k = "{k}",
		r = "{r}",
		hRat = "{hr}",
		d = "{d}",
		i = "{i}"
	output:
		homs = temp("../simulations/homologies/homologies_{genome}_{desc}_k{k}_hr{hr}_c{c}_u{u}_de{d}_in{i}" + \
			"_rep{r}.txt"),
		bench = "../benchmarks/benchFindThoms_{genome}_{desc}_k{k}_hr{hr}_c{c}_u{u}_de{d}_in{i}_rep{r}.txt"
	wildcard_constraints:
		genome = "\w+",
	shell:
		"/usr/bin/time -v src/FindThoms -p {input.rds} -s {input.txt} -k {params.k} -r {params.hRat} -c {params.c} -u " + \
		"{params.u} -d {params.d} -i {params.i} -N > {output.homs} 2> {output.bench}"

rule searchHumanHomologies:
	input:
		pttn = "../simulations/reads/humanChr20_cR103_ep{errorPattern}_s{seed}_rr{i}.fasta",
		txt = "../simulations/genomes/{refName}.fasta"
	params:
		c = "{c}",
		u = "{u}",
		thres = "{t}",
		k = "{k}"
	output:
		homs = "../simulations/homologies/homologies_{refName}_ep{errorPattern}_s{seed}_rr{i}_k{k}_c{c}_u{u}_t{t}.txt",
		bench = "../benchmarks/benchFindThoms_{refName}_ep{errorPattern}_s{seed}_rr{i}_k{k}_c{c}_u{u}_t{t}.txt"
	shell:
		"/usr/bin/time -v src/FindThoms -p {input.pttn} -s {input.txt} -k {params.k} -c {params.c} -u {params.u} -t " + \
		"{params.thres} > {output.homs} 2> {output.bench}"

rule searchHomologies:
	input:
		pttn = "../simulations/reads/reads_tl{tl}_rid{ri}_ch{chem}_ep{erR}_{i}.fasta",
		txt = "../simulations/text_rl{rl}_rid{ri}_m{desc}.fasta"
	params:
		c = "{c}",
		u = "{u}",
		thres = "{t}",
		k = "{k}"
	output:
		homs = "../simulations/homologies/homologies_rl{rl}_rid{ri}_m{desc}_tl{tl}_ch{chem}_ep{erR}_ri{i}_k{k}_c{c}_u{u}_t{t}.txt",
		bench = "../benchmarks/benchFindThoms_rl{rl}_rid{ri}_m{desc}_tl{tl}_ch{chem}_ep{erR}_ri{i}_k{k}_c{c}_u{u}_t{t}.txt"
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
