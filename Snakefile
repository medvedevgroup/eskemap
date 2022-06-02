configfile: 'config.yaml'

rule all:
	input:
		expand("simulations/simSeqs_s{s}_n{n}_l{l}_m{m}_i{i}_d{d}.ssv", s=config['randomSeed'], n=config['nbSimSeqs'], l=config['simSeqLen'], m=\
			config['mutationRates'])

rule simSeqs:
	output:
		"simSeqs_s{s}_n{n}_l{l}_m{m}_i{i}_d{d}.ssv"
	shell:
		"python3 SimSeqPairs.py -s {wildcards.s} -n {wildcards.n} -l {wildcards.l} -m {wildcards.m} -i {wildcards.i} -d {wildcards.d} > {output}"
