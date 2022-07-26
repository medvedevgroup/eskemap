#ifndef IO_HPP
#define IO_HPP

#include "Measures.h"
#include "Thomology.h"

#define OPTIONS "a:b:ilh"
#define T_HOM_OPTIONS "p:s:k:r:c:u:t:nh"
#define MIN_PARAM_NB 6
#define MAX_RATIO 1.0
#define NORM_FLAG_DEFAULT false

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
	cerr << "FindThoms [-hn] [-p PATTERN] [-s TEXT] [-k KMER_LEN] [-r HASH_RATIO] [-c COM_HASH_WGHT] [-u UNI_HASH_WGHT] [-t HOM_THRES]" << endl << endl;
	cerr << "Find sketch-based pattern homology in text." << endl << endl;
	cerr << "Required parameters:" << endl;
	cerr << "   -p   --pattern  Pattern sequence" << endl;
	cerr << "   -s   --text     Text sequence" << endl << endl;
	cerr << "Optional parameters with required argument:" << endl;
	cerr << "   -k   --ksize             K-mer length to be used for sketches" << endl;
	cerr << "   -r   --hashratio         FracMin hash ratio to be used for sketches" << endl;
	cerr << "   -c   --commonhashweight  Weight to reward common hashes (default " << DEFAULT_WEIGHT << ")" << endl;
	cerr << "   -u   --uniquehashweight  Weight to punish unique hashes (default " << DEFAULT_WEIGHT << ")" << endl;
	cerr << "   -t   --hom_thres         Homology threshold (default " << T << ")" << endl;
	cerr << "Optional parameters without argument:" << endl;
	cerr << "   -n   --normalize  Normalize scores by length" << endl;
	cerr << "   -h   --help       Display this help message" << endl;
}

//This function parses the program parameters. Returns false if given arguments are not valid
const bool prsArgs(int& nArgs, char** argList, string& seqa, string& seqb, Measure& msr);

//This function parses the program parameters. Returns false if given arguments are not valid
const bool prsArgs(int& nArgs, char** argList, string& pat, string& txt, uint32_t& k, double& hFrac, uint32_t& cw, uint32_t& uw, int32_t& tThres, 
	bool& norm);

//This function reads a file in FASTA format and returns true on success
const bool readFASTA(const string& filePath, string& seq);

#endif