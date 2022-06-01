#ifndef IO_HPP
#define IO_HPP

#include "Measures.h"

#define OPTIONS "a:b:ih"
#define MIN_PARAM_NB 6

//This function prints usage infos
inline void dspHlp(){
	cerr << "CalcSim [-hi] [-a SEQ1] [-b SEQ2]" << endl << endl;
	cerr << "Calculating sequence similarity based on their sketches." << endl << endl;
	cerr << "Required parameters:" << endl;
	cerr << "   -a   --seqa  First input sequence" << endl;
	cerr << "   -b   --seqb  Second input sequence" << endl << endl;
	cerr << "Required parameters without argument:" << endl;
	cerr << "   -i   --intersim  Calculate similarity based on intersection measure" << endl;
	cerr << "Optional parameters without argument:" << endl;
	cerr << "   -h   --help      Display this help message" << endl;
}

//This function parses the program parameters. Returns false if given arguments are not valid
const bool prsArgs(int& nArgs, char** argList, string& seqa, string& seqb, Measure& msr);

#endif