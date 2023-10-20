#include <string>
#include <iostream>
#include <cstdint>
#include <unordered_map>

#include "../src/Sketch.cpp"
#include "../src/IO.cpp"

#define INDEX_DEFAULT_DUMP_FILE "indexDump.idx"

unordered_map<uint64_t, char> bLstmers;

int main(int argc, char **argv){
	//The number of occurrences of a k-mer
	int nbOcc;
	//Maximum number of times a k-mer is allowed to occur
	int maxAb;
	//Input file name
	string seqFile;
	//Input sequence
	string seq;
	//An index option struct
	mm_idxopt_t iopt;
	//A mapping options struct
	mm_mapopt_t mopt;
	//An index reader
	mm_idx_reader_t *r;
	//A pointer to the sequence's index
	const mm_idx_t *idx;
	//Input sequence's sketch
	Sketch seqSketch;
	//Hash table to save which k-mer hashes we have seen already
	unordered_map<uint64_t, char> seenHashes;

	//Set index options to default
	mm_set_opt(0, &iopt, &mopt);

	//Parse arguments
	if(argc < 5){
		cerr << "Usage: CreateKmerBlacklist [SeqFile] [kSize] [wSize] [maxKmerAbundance]" << endl;
		return 1;
	}

	seqFile = argv[1];
	iopt.k = atoi(argv[2]);
	iopt.w = atoi(argv[3]);
	maxAb = atoi(argv[4]);
	//Open an index reader
	r = mm_idx_reader_open(seqFile.c_str(), &iopt, INDEX_DEFAULT_DUMP_FILE);

	//Check if reader could be opened successfully
	if(r == NULL){
		cerr << "ERROR: Text sequence file could not be read" << endl;
		return 1;
	}

	//Construct index 
	if((idx = mm_idx_reader_read(r, 1)) == 0){
		cerr << "ERROR: Text index cannot be read" << endl;
		return 1;
	}

	//Read sequence file
	readFASTA(seqFile, seq);
	//Build sketch
	seqSketch = buildMiniSketch(seq, idx->k, idx->w, bLstmers);

	//Iterate over k-mer hashes in sketch
	for(Sketch::const_iterator i = seqSketch.begin(); i != seqSketch.end(); ++i){
		if(!seenHashes.contains(*i)){
			seenHashes[*i] = 1;
			const uint64_t *idx_p = mm_idx_get(idx, *i, &nbOcc);
			if(nbOcc > maxAb)	cout << *i << endl;
		}
	}

	return 0;
}
