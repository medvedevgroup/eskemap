# ESKEMAP - Exact sketch-based read mapping

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
FindThoms [-hn] [-p PATTERN_FILE] [-s TEXT_FILE] [-k KMER_LEN] [-r HASH_RATIO] [-b BLACKLIST] [-c COM_HASH_WGHT] [-u UNI_HASH_WGHT] [-t HOM_THRES] [-d DECENT] [-i INTERCEPT] [-N NESTING]

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

**TODO...**

### Data

The T2T reference assembly of human chromosome Y (Accession number NC_060948.1) was downloaded from [NCBI](https://www.ncbi.nlm.nih.gov).

### Reproduction Workflow

**TODO...**

### Notebook Analysis

**TODO...**

## Support

For any questions, feedback or problem, please feel free to file an issue or contact the developers via email and we will get back to you as soon as possible.

## Licenses

* minimap2 is licensed under the MIT license (https://github.com/lh3/minimap2)

* ESKEMAP is GNU GPLv3 licensed 
