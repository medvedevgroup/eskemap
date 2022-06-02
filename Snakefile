configfile: 'config.yaml'

rule all:
	input:
		expand("../simulations/scores_s{s}_n{n}_l{l}_m{m}_i{m}_d{m}.txt", s=config['randomSeed'], n=config['nbSimSeqs'], l=config['simSeqLen'], m=\
			config['mutationRates'])

rule simSeqs:
	output:
		"../simulations/simSeqs_s{s}_n{n}_l{l}_m{m}_i{i}_d{d}.ssv"
	shell:
		"python3 SimSeqPairs.py -s {wildcards.s} -n {wildcards.n} -l {wildcards.l} -m {wildcards.m} -i {wildcards.i} -d {wildcards.d} > {output}"

rule measureSimilarity:
	input:
		"../simulations/simSeqs_{simDesc}.ssv"
	output:
		"../simulations/scores_{simDesc}.txt"
	run:
		for p in open(input[0], 'r'):
			seqs = p.strip().split(' ')
			shell("src/CalcSim -i -a {seqs[0]} -b {seqs[1]} >> {output}")
