configfile: 'config.yaml'

from random import randrange, seed
from sys import maxsize
from glob import glob

SUB_ERR = 0.002 * 6 / 110.
INS_ERR = 0.002 * 50 / 110.
DEL_ERR = 0.002 * 54 / 110.

#Initialize random number generator with global seed for this workflow
seed(config['globalSeed'])

def enumerateEdlibRes(wcs):
	smallestId, largestId = [int(n) for n in wcs.rng.split('-')]

	return expand("simulations/edlibMappings/{d}_ri{i}.er", d=wcs.desc, i=[i for i in range(smallestId, largestId + 1)])

READ_SEED = randrange(maxsize)

rule all:
	input:
		"simulations/edlibMappings/%s_sr%.19f_dr%.19f_i%.19f_sd%d_lmn100_lmx1000000_lavg9000_ls7000_dp10_ri0-69400.er" \
		%(config['ref'], SUB_ERR, DEL_ERR, INS_ERR, READ_SEED),
		"simulations/homologies/homologies_%s_sr%.19f_dr%.19f_i%.19f_sd" %(config['ref'], SUB_ERR, DEL_ERR, INS_ERR) + \
		"%d_lmn100_lmx1000000_lavg9000_ls7000_dp10_rm20_k15_w10_c1_u1_de%.8f_in%.13f.txt" %(READ_SEED, config['eskemapDecent'], \
			config['eskemapIntercept']),
		"benchmarks/benchEskemap_%s_sr%.19f_dr%.19f_i%.19f_sd%d_lmn100_"  %(config['ref'], SUB_ERR, DEL_ERR, INS_ERR, READ_SEED) + \
		"lmx1000000_lavg9000_ls7000_dp10_rm20_k15_w10_c1_u1_de%.8f_in%.13f_rep0.txt" %(config['eskemapDecent'], config\
			['eskemapIntercept']),
		"simulations/minimap2Res/%s_sr%.19f_dr%.19f_i%.19f_sd%d_lmn100_lmx1000000_lavg9000_ls7000_dp10_rm20_k15.paf.gz" \
		%(config['ref'], SUB_ERR, DEL_ERR, INS_ERR, READ_SEED),
		"benchmarks/benchMinimap2ApprxMppng_%s_sr%.19f_dr%.19f_i%.19f_sd" %(config['ref'], SUB_ERR, DEL_ERR, INS_ERR) + \
		"%d_lmn100_lmx1000000_lavg9000_ls7000_dp10_rm20_k15_rep0.txt" %READ_SEED,
		"simulations/Winnowmap2Res/%s_sr%.19f_dr%.19f_i%.19f_sd%d_lmn100_lmx1000000_lavg9000_ls7000_dp10_rm20_k15.paf.gz" \
		%(config['ref'], SUB_ERR, DEL_ERR, INS_ERR, READ_SEED),
		"benchmarks/benchWinnowmap2ApprxMppng_%s_sr%.19f_dr%.19f_i%.19f_sd" %(config['ref'], SUB_ERR, DEL_ERR, INS_ERR) + \
		"%d_lmn100_lmx1000000_lavg9000_ls7000_dp10_rm20_k15_rep0.txt" %READ_SEED,
		"simulations/blastRes/subs_Edlib_%s_sr%.19f_dr%.19f_i%.19f_sd%d" %(config['ref'], SUB_ERR, DEL_ERR, INS_ERR, READ_SEED) + \
			"_lmn100_lmx1000000_lavg9000_ls7000_dp10_ri0-69400_rm20_e0.01.tsv",
		"simulations/blastRes/subs_ESKEMAP_%s_sr%.19f_dr%.19f_i%.19f_sd" %(config['ref'], SUB_ERR, DEL_ERR, INS_ERR) + \
			"%d_lmn100_lmx1000000_lavg9000_ls7000_dp10_rm20_k15_w10_c1_u1_de%.8f_in%.13f_e0.01.tsv" %(READ_SEED, config[\
				'eskemapDecent'], config['eskemapIntercept']),
		expand("simulations/blastRes/subs_{t}_%s_sr%.19f_dr%.19f_i%.19f_sd" %(config['ref'], SUB_ERR, DEL_ERR, INS_ERR) + \
			"%d_lmn100_lmx1000000_lavg9000_ls7000_dp10_rm20_k15_e0.01.tsv" %READ_SEED, t=config['compTools'])

rule blastPairwiseMultFasta:
	input:
		"simulations/mappedAreas/{desc}.fasta"
	output:
		"simulations/blastRes/{desc}_e10.tsv"
	shell:
		"python3 scripts/BlastPairwiseMultiFasta.py -f {input} -o {output}"

rule filterBlastRes:
	input:
		"simulations/blastRes/{desc}_e10.tsv"
	params:
		"{ev}"
	output:
		"simulations/blastRes/{desc}_e{ev}.tsv"
	wildcard_constraints:
		ev = "0\.[0-9]*"#"[0-9]?\.[0-9]*"
	run:
		ofile = open(output[0], 'w')

		for l in open(input[0], 'r'):
			if l.startswith("Results"):
				ofile.write(l)
			elif float(l.split('\t')[6]) <= float(params[0]):
				ofile.write(l)

		ofile.close()

rule getCompResRefSubstrings:
	input:
		ref = "simulations/genomes/{genome}.fasta",
		res = "simulations/{toolPrefix}map2Res/{genome}_sr{desc}.paf.gz"
	output:
		"simulations/mappedAreas/subs_{toolPrefix}map2_{genome}_sr{desc}.fasta"
	shell:
		"python3 scripts/getRefSubstringsFromPAFres.py -r {input.ref} -p {input.res} -o {output}"

rule getESKEMAPresRefSubstrings:
	input:
		ref = "simulations/genomes/{genome}.fasta",
		bl = "%s.txt" %config['kmerBlacklistName'],
		res = "simulations/homologies/homologies_{genome}_sr{desc}_k{k}_w{w}_c{mdesc}.txt"
	params:
		k = "{k}",
		w = "{w}"
	output:
		"simulations/mappedAreas/subs_ESKEMAP_{genome}_sr{desc}_k{k}_w{w}_c{mdesc}.fasta"
	shell:
		"python3 scripts/getRefSubstringsFromESKEMAPres.py -r {input.ref} -b {input.bl} -k {params.k} -w {params.w} -e " + \
		"{input.res} -o {output}"

rule getEdlibResRefSubstrings:
	input:
		ref = "simulations/genomes/{genome}.fasta",
		res = "simulations/edlibMappings/{genome}_sr{desc}.er"
	params:
		"{thres}"
	output:
		"simulations/mappedAreas/subs_Edlib_{genome}_sr{desc}_rm{thres}.fasta"
	shell:
		"python3 scripts/getRefSubstringsFromEdlibRes.py -r {input.ref} -m {input.res} -t {params} -o {output}"

rule saveWinnowmap2Result:
	input:
		"simulations/Winnowmap2Res/{genome}_sr{desc}_k{k}_rep0.{frmt}.gz"
	output:
		"simulations/Winnowmap2Res/{genome}_sr{desc}_k{k}.{frmt}.gz"
	wildcard_constraints:
		k = "[0-9]+"
	shell:
		"mv {input} {output}"

rule runApprxMppngWinnowmap2onRealGenomeFASTA:
	input:
		ref = "simulations/genomes/{genome}.fasta",
		qry = "simulations/reads/{genome}_sr{desc}.fasta",
		cnts = "simulations/repKmers_k{k}_{genome}.txt"
	params:
		k = "{k}",
		r = "{r}"
	output:
		res = temp("simulations/Winnowmap2Res/{genome}_sr{desc}_k{k}_rep{r}.paf.gz"),
		bench = "benchmarks/benchWinnowmap2ApprxMppng_{genome}_sr{desc}_k{k}_rep{r}.txt"
	wildcard_constraints:
		r = "[0-9]+"
	shell:
		"/usr/bin/time -v %s -W {input.cnts} -k {params.k} {input.ref} {input.qry} " %config['WinnowmapBin'] + \
		"2> {output.bench} | gzip -3 > {output.res}"

rule printCounts:
	input:
		"simulations/merylDB_k{k}_{desc}"
	output:
		temp("simulations/repKmers_k{k}_{desc}.txt")
	shell:
		"%s print greater-than distinct=0.9998 {input} > {output}" %config['merylBin']

rule countGenomeKmers:
	input:
		"simulations/genomes/{genome}.fasta"
	params:
		"{k}"
	output:
		temp(directory("simulations/merylDB_k{k}_{genome}"))
	shell:
		"%s count k={params} output {output} {input}" %config['merylBin']

rule saveMinimap2Result:
	input:
		"simulations/minimap2Res/{genome}_sr{desc}_k{k}_rep0.{frmt}.gz"
	output:
		"simulations/minimap2Res/{genome}_sr{desc}_k{k}.{frmt}.gz"
	wildcard_constraints:
		k = "[0-9]+"
	shell:
		"mv {input} {output}"

rule runApprxMppngMinimap2onRealGenomePacBioFASTA:
	input:
		ref = "simulations/genomes/{genome}.fasta",
		qry = "simulations/reads/{genome}_sr{desc}.fasta"
	params:
		k = "{k}",
		r = "{r}"
	output:
		res = temp("simulations/minimap2Res/{genome}_sr{desc}_k{k}_rep{r}.paf.gz"),
		bench = "benchmarks/benchMinimap2ApprxMppng_{genome}_sr{desc}_k{k}_rep{r}.txt"
	wildcard_constraints:
		r = "[0-9]+"
	shell:
		"/usr/bin/time -v %s {input.ref} {input.qry} -k {params.k} 2> {output.bench} | gzip -3 > {output.res}" \
		%config['minimap2Bin']

rule saveFindThomsResult:
	input:
		"simulations/homologies/homologies_{genome}_{desc}_k{k}_{smp}_c{c}_u{u}_de{d}_in{i}_rep0.txt"
	output:
		"simulations/homologies/homologies_{genome}_{desc}_k{k}_{smp}_c{c}_u{u}_de{d}_in{i}.txt"
	wildcard_constraints:
		i = "-?[0-9]+\.?[0-9]*"
	shell:
		"mv {input} {output}"

rule filterReads:
	input:
		e = "simulations/edlibMappings/{desc}_ri0-69400.er",
		r = "simulations/reads/{desc}.fasta"
	params:
		"{rm}"
	output:
		"simulations/reads/{desc}_rm{rm}.fasta"
	shell:
		"python3 scripts/FilterReads.py -e {input.e} -r {input.r} -m {params} -o {output}"

rule searchMinimapSketchReadHomologies:
	input:
		rds = "simulations/reads/{genome}_{desc}.fasta",
		txt = "simulations/genomes/{genome}.fasta",
		bl = "%s.txt" %config['kmerBlacklistName']
	params:
		c = "{c}",
		u = "{u}",
		k = "{k}",
		r = "{r}",
		w = "{w}",
		d = "{d}",
		i = "{i}"
	output:
		homs = temp("simulations/homologies/homologies_{genome}_{desc}_k{k}_w{w}_c{c}_u{u}_de{d}_in{i}" + \
			"_rep{r}.txt"),
		bench = "benchmarks/benchEskemap_{genome}_{desc}_k{k}_w{w}_c{c}_u{u}_de{d}_in{i}_rep{r}.txt"
	wildcard_constraints:
		genome = "\w+",
	shell:
		"/usr/bin/time -v src/eskemap -p {input.rds} -s {input.txt} -k {params.k} -c {params.c} -u " + \
		"{params.u} -d {params.d} -i {params.i} -w {params.w} -b {input.bl} -N > {output.homs} 2> {output.bench}"

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
	wildcard_constraints:
		rdId = "[0-9]+"
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
	run:
		ofile = open(output[0], 'w')

		for f in input:
			ofile.write(f + '\n')
			
			for l in open(f, 'r'):
				ofile.write(l)

		ofile.close()
