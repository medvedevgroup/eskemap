configfile: 'config.yaml'

from random import randrange, seed
from sys import maxsize

SUB_ERR = 0.002 * 6 / 110.
INS_ERR = 0.002 * 50 / 110.
DEL_ERR = 0.002 * 54 / 110.

#Initialize random number generator with global seed for this workflow
seed(config['globalSeed'])

def enumerateEdlibRes(wcs):
	smallestId, largestId = [int(n) for n in wcs.rng.split('-')]

	return expand("../simulations/edlibMappings/{d}_ri{i}.er", d=wcs.desc, i=[i for i in range(smallestId, largestId + 1)])

READ_SEED = randrange(maxsize)

rule all:
	"simulations/edlibMappings/t2thumanChrY_sr%.19f_dr%.19f_i%.19f_sd%d_lmn100_lmx1000000_lavg9000_ls7000_dp10_ri0-69401.er" \
	%(SUB_ERR, DEL_ERR, INS_ERR, READ_SEED),


#Simulate read set
rule simReadsOwnScript:
	input:
		ref = "simulations/genomes/{genome}.fasta"
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
		rds = "simulations/reads/{genome}_sr{subR}_dr{delR}_i{insR}_sd{sd}_lmn{lMin}_lmx{lMax}_lavg{lMn}_ls{lStd}_dp{dp}.fasta"
	wildcard_constraints:
		dp = "[0-9]+"
	shell:
		"python3 scripts/simReads.py -dp {params.dp} -lmn {params.lMin} -lmx {params.lMax} -lavg {params.lMean} -ls " + \
		"{params.lStd} -r {input} -sr {params.subR} -dr {params.delR} -ir {params.insR} -sd {params.sd} -o {output}"

rule divideReads:
	input:
		"simulations/reads/{rdFileName}.fasta"
	params:
		"{rdId}"
	output:
		"simulations/reads/{rdFileName}_ri{rdId}.fasta"
	shell:
		"python3 scripts/getSeq.py -s {input} -i {params} -o {output}"

rule runEdlib:
	input:
		ref = "simulations/genomes/{genome}.fasta",
		qry = "simulations/reads/{genome}_sr{rdDesc}_ri{rdId}.fasta"
	output:
		temp("simulations/edlibMappings/{genome}_sr{rdDesc}_ri{rdId, [0-9]+}.er")
	shell:
		"FindSimSeqs/FindSimSeqs {input.qry} {input.ref} > {output}"

rule mergeEdlibRes:
	input:
		enumerateEdlibRes
	output:
		"simulations/edlibMappings/{desc}_ri{rng, [0-9]+-[0-9]+}.er"
	shell:
		"for f in {input}; do echo $f; cat $f; done > {output}"
