# ESKEMAP - Exact SKEtch-based read MAPping

This directory contains the source code of ESKEMAP and the documentation of experiments we performed in our paper.

## Description

ESKEMAP is an algorithm for sketch-based, all-hits read mapping of genomic sequences. Using the sketches of a read and a reference genome, it finds all positions inside the reference genome's sketch that score high enough with respect to a similarity measure and a threshold depending on the read's length.

## Installation

ESKEMAP can be used with the sketching approach implemented in [minimap2](https://github.com/lh3/minimap2). A modified version of its source code is located inside the subdirectory *minimap2*. It can be compiled using cmake.

```
cd <corer_directory>/minimap2
make
```

Afterwards, ESKEMAP can be compiled from the directory *src*.

```
cd <corer_directory>/src
make
```

## Usage

```
./eskemap
```

displays the command line interface:
```
eskemap [-hn] [-p PATTERN_FILE] [-s TEXT_FILE] [-k KMER_LEN] [-r HASH_RATIO] [-b BLACKLIST] [-c COM_HASH_WGHT] [-u UNI_HASH_WGHT] [-t HOM_THRES] [-d DECENT] [-i INTERCEPT] [-N NESTING]

Find sketch-based pattern similarity in text.

Required parameters:
   -p   --pattern  Pattern sequences file (FASTA format)
   -s   --text     Text sequence file (FASTA format)

Optional parameters with required argument:
   -k   --ksize             K-mer length to be used for sketches (default 9)
   -w   --windowsize        Window size for minimizer sketching approach (default 10)
   -r   --hashratio         FracMin hash ratio to be used for sketches (default 0.1)
   -b   --blacklist         File containing hashes to ignore for sketch calculation
   -c   --commonhashweight  Weight to reward common hashes (default 1)
   -u   --uniquehashweight  Weight to punish unique hashes (default 1)
   -t   --hom_thres         Homology threshold (default 0)
   -d   --decent            Decent required for dynamic threshold selection
   -i   --intercept         Intercept required for dynamic threshold selection
Optional parameters without argument:
   -n   --normalize  Normalize scores by length
   -N   --nesting    Output nested homologies
   -h   --help       Display this help message
```

A common program call of ESKEMAP would look like:

```
eskemap -p reads.fa -s reference.fa -k 15 -w 10 -b highFreqKmers.txt -d 0.8 -i 12 -N
```

Here, we let ESKEMAP find mappings for all read sequences stored inside *reads.fa* (in FASTA format) to the sequence stored inside *reference.fa* (in FASTA format). Minimizer sketches (see [Roberts *et al.*](https://doi.org/10.1093/bioinformatics/bth408) and [Schleimer *et al.*](https://doi.org/10.1145/872757.872770) for details) are used with *k*=15 and window size *w*=10. All *k*-mers listed in *highFreqKmers.txt* are excluded from all sketches. Score thresholds are calculated with respect to a read's length using a line function with decent 0.8 and intercept 12. All mappings including nested ones (`-N`) are reported.

## Experiments

This section documents how the experiments described in our paper can be reproduced.

### Requirements

Our experiments are partly documented as a [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow that allows to rerun those parts using exact program calls.
The remaining steps of our experiment including postprocessing for result analysis and plot generation can be deduced/rerun from a [Jupyter Notebook](https://jupyter.org). For a reproduction of the whole experiment, both tools need to be installed on your local system. Additionally, the Python package [Biopython](https://biopython.org) and a running version of [minimap2](https://github.com/lh3/minimap2/tree/master) and [Winnowmap2](https://github.com/marbl/Winnowmap) are also required.

We also used the API of [Edlib](https://github.com/Martinsos/edlib#api-documentation) and [parasail](https://github.com/jeffdaily/parasail) to implement a small script to calculate read mapping positions on the basis of alignments. After Edlib and parasail are installed on your system, the script can be installed by changing into the subdirectory *FindSimSeqs* and executing `make`.

```
cd FindSimSeqs
make
```

### Data

#### Human Chromosome Y

The T2T reference assembly of human chromosome Y (Accession number NC_060948.1) was downloaded from [NCBI](https://www.ncbi.nlm.nih.gov).

#### Simulated Reads

In order to create the exact same set of reads we used for our experiments, run:

```
mkdir -p simulations/reads
python3 scripts/simReads.py -dp 10 -lmn 100 -lmx 1000000 -lavg 9000 -ls 7000 -r simulations/genomes/t2thumanChrY.fasta -sr 0.00010909090909090909 -dr 0.0009818181818181818 -ir 0.0009090909090909091 -sd 7361077429744071834 -o simulations/reads/t2thumanChrY_sr0.00010909090909090909_dr0.0009818181818181818_i0.0009090909090909091_sd7361077429744071834_lmn100_lmx1000000_lavg9000_ls7000_dp10.fasta
```

This will create a subdirectory structure *simulations/reads* containing the unfiltered, whole set of reads stored as a FASTA file. 

For filtering the read sets according to edlib's mapping results, these results first have to be generated. This can be done using the snakemake workflow:

```
snakemake simulations/edlibMappings/t2thumanChrY_sr0.00010909090909090909_dr0.0009818181818181818_i0.0009090909090909091_sd7361077429744071834_lmn100_lmx1000000_lavg9000_ls7000_dp10_ri0-69400.er
```

How to create a reads file only containing those reads for which edlib could find only up to 20 different, non-overlapping mapping positions is documented inside the [Jupyter Notebook](https://jupyter.org) *Experiments.ipynb* (Section *Read Filtering*). Rerunning the respective notebook part will create a file ending with ...`_dp10_rm20.fasta` inside the directory *simulations/reads*.

### Reproduction Workflow

Exact program calls of each program we used for mapping the reads are documented inside a [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow. The workflow may even be run to reproduce the whole experiment if all dependencies are satisfied (see [Requirements](#Requirements)) and the necessary input data is provided (see [Data](#Data)).

Before running the workflow a few configuration steps need to be done to ensure that program binaries and input files are found.

**TODO...**

Afterwards, the workflow can be run by just typing:

```
snakemake
```

**TODO...**

### Notebook Analysis

**TODO...**

## Support

For any questions, feedback or problem, please feel free to file an issue or contact the developers via email and we will get back to you as soon as possible.

## Licenses

* minimap2 is licensed under the MIT license (https://github.com/lh3/minimap2)

* edlib is licensed under the MIT license (https://github.com/Martinsos/edlib)

* parasail is licensed by the Battelle Memorial Institute

* ESKEMAP is GNU GPLv3 licensed 
