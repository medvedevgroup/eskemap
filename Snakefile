configfile: 'config.yaml'

NB_RANDSEQS=(config['nbCpys'] + 1) * config['nbSimSeqs']

rule all:
	input:
		expand("../simulations/homologies_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{m}_d{m}.txt", gn=config['nbSimSeqs'], rn=NB_RANDSEQS, gl=\
			config['geneLen'], rl=config['randSeqLen'], o=config['nbCpys'], m=config['mutationRates'] + [0])#,
		#expand("../simulations/scores_mes{mes}_n{n}_l{l}_m{m}_i{m}_d{m}.txt", mes=config['simMeasure'], n=\
		#	config['nbSimSeqs'], l=config['simSeqLen'], m=config['mutationRates'])

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
		"../simulations/simGeneSeqs_n{n}_l{l}_o{o}_m{m}_i{i}_d{d}.ssv"
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
		"../simulations/simGeneSeqs_n{gn}_l{gl}_o{o}_m{m}_i{i}_d{d}.ssv",
		"../simulations/randSeqs_n{rn}_l{rl}_m{m}_i{i}_d{d}.ssv"
	output:
		"../simulations/searchPairs_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{i}_d{d}.txt"
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
		"../simulations/searchPairs_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{i}_d{d}.txt"
	output:
		"../simulations/homologies_gn{gn}_rn{rn}_gl{gl}_rl{rl}_o{o}_m{m}_i{i}_d{d}.txt"
	run:
		i = 0
		for l in open(input[0], 'r'):
			i += 1
			shell("echo Pair {i} >> {output}")
			l = l.strip()
			pattern, text = l.split(' ')
			shell("python3 scripts/FindThoms.py -p {pattern} -s {text} >> {output}")
