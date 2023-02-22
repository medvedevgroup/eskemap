#ifndef IO_HPP
#define IO_HPP

#include "Measures.h"
#include "Thomology.h"

#define OPTIONS "a:b:ilh"
#define T_HOM_OPTIONS "p:s:k:r:b:c:u:t:d:i:nNh"
#define MIN_PARAM_NB 6
#define MAX_RATIO 1.0
#define NORM_FLAG_DEFAULT false
#define NESTING_FLAG_DEFAULT true
#define PATTERN_BATCH_SIZE 250000
#define STRING_BUFFER_SIZE_DEFAULT 50

//This function prints usage infos
inline void dspHlp(){
	cerr << "CalcSim [-hil] [-a SEQ1] [-b SEQ2]" << endl << endl;
	cerr << "Calculating sequence similarity based on their sketches." << endl << endl;
	cerr << "Required parameters:" << endl;
	cerr << "   -a   --seqa  First input sequence" << endl;
	cerr << "   -b   --seqb  Second input sequence" << endl << endl;
	cerr << "Required parameters without argument:" << endl;
	cerr << "   -i   --intersim     Calculate similarity based on intersection measure" << endl;
	cerr << "   -l   --algnHshsSim  Calculate similarity based on aligning hashes measure" << endl;
	cerr << "Optional parameters without argument:" << endl;
	cerr << "   -h   --help      Display this help message" << endl;
}

//This function prints usage infos
inline void dsHlp(){
	cerr << "FindThoms [-hn] [-p PATTERN_FILE] [-s TEXT_FILE] [-k KMER_LEN] [-r HASH_RATIO] [-b BLACKLIST] [-c COM_HASH_WGHT] \
[-u UNI_HASH_WGHT] [-t HOM_THRES] [-d DECENT] [-i INTERCEPT] [-N NESTING]" << endl << endl;
	cerr << "Find sketch-based pattern homology in text." << endl << endl;
	cerr << "Required parameters:" << endl;
	cerr << "   -p   --pattern  Pattern sequences file (FASTA format)" << endl;
	cerr << "   -s   --text     Text sequence file (FASTA format)" << endl << endl;
	cerr << "Optional parameters with required argument:" << endl;
	cerr << "   -k   --ksize             K-mer length to be used for sketches (default " << K << ")" << endl;
	cerr << "   -r   --hashratio         FracMin hash ratio to be used for sketches (default " << HASH_RATIO << ")" << endl;
	cerr << "   -b   --blacklist         File containing hashes to ignore for sketch calculation" << endl;
	cerr << "   -c   --commonhashweight  Weight to reward common hashes (default " << DEFAULT_WEIGHT << ")" << endl;
	cerr << "   -u   --uniquehashweight  Weight to punish unique hashes (default " << DEFAULT_WEIGHT << ")" << endl;
	cerr << "   -t   --hom_thres         Homology threshold (default " << T << ")" << endl;
	cerr << "   -d   --decent            Decent required for dynamic threshold selection" << endl;
	cerr << "   -i   --intercept         Intercept required for dynamic threshold selection" << endl;
	cerr << "Optional parameters without argument:" << endl;
	cerr << "   -n   --normalize  Normalize scores by length" << endl;
	cerr << "   -N   --nesting    Output nested homologies" << endl;
	cerr << "   -h   --help       Display this help message" << endl;
}

//This function outputs all given t-homologies
inline void outputHoms(const vector<Thomology>& homs, const bool& norm, const uint32_t& pLen){
	//Iterate over t-homologies
	for(vector<Thomology>::const_iterator h = homs.begin(); h != homs.end(); ++h){
		//Normalize score before reporting if requested
		if(norm){
			cout << "i: " << get<0>(*h) << " j: " << get<1>(*h) << " score: " << (double) get<2>(*h) / max(pLen, get<1>(*h) - 
				get<0>(*h) + 1) << endl;
		} else{
			cout << "i: " << get<0>(*h) << " j: " << get<1>(*h) << " score: " << get<2>(*h) << endl;
		}
	}
}

//This function parses the program parameters. Returns false if given arguments are not valid
const bool prsArgs(int& nArgs, char** argList, string& seqa, string& seqb, Measure& msr);

//This function parses the program parameters. Returns false if given arguments are not valid
const bool prsArgs(int& nArgs, char** argList, string& pFl, string& tFl, uint32_t& k, double& hFrac, string& blFl, uint32_t& cw, 
	float& uw, float& tThres, bool& norm, float& dec, float& inter, bool& noNesting);

//This function reads a file in FASTA format and returns true on success
const bool readFASTA(const string& filePath, string& seq);

//This function reads in batches of FASTA sequence entries from file and transforms them into sketches. Returns false if end of file
//was reached.
const bool lPttnSks(ifstream& fStr, const uint32_t& k, const double& hFrac, const unordered_map<uint64_t, char>& bLstmers, 
	vector<tuple<string, uint32_t, Sketch>>& pSks);

//This function reads in batches of FASTA sequence entries from file and transforms them into minimap sketches. Returns false if end
// of file was reached.
const bool lMiniPttnSks(ifstream& fStr, const uint32_t& k, const uint32_t& w, const unordered_map<uint64_t, char>& blmers, 
	vector<tuple<string, uint32_t, Sketch>>& pSks);

//This function reads 64-bit numbers from file and returns them as a hash table
const unordered_map<uint64_t, char> readBlstKmers(const string& fname);

#endif