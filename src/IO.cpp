#include <getopt.h>

#include "IO.h"

//This function parses the program parameters. Returns false if given arguments are not valid
const bool prsArgs(int& nArgs, char** argList, string& seqa, string& seqb, Measure& msr){
	bool seqaGvn = false, seqbGvn = false;
	int option_index = 0, a;

	//Check if enough arguments are given at all
	if(nArgs < MIN_PARAM_NB) return false;

	static struct option long_options[] = {
        {"seqa",        required_argument,  0, 'a'},
        {"seqb",        required_argument,  0, 'b'},
        {"intersim",    no_argument,        0, 'i'},
        {"algnHshsSim", no_argument,        0, 'l'},
        {"help",        no_argument,        0, 'h'},
        {0,             0,                  0,  0 }
    };

    //Parse all parameters given
	while ((a = getopt_long(nArgs, argList, OPTIONS, long_options, &option_index)) != -1){
		//Assign parameter values
		switch(a){
			case 'a':
				//Save input sequence
				seqa = optarg;
				//Note that we have seen one input sequence
				seqaGvn = true;
				break;
			case 'b':
				//Save input sequence
				seqb = optarg;
				//Note that we have seen one input sequence
				seqbGvn = true;
				break;
			case 'i':
				//Check if another measure was also given as a parameter
				if(msr != none && msr != intersec){
					cerr << "ERROR: Only one measure can be calculated at a time!" << endl;

					return false;
				}

				msr = intersec;
				break;
			case 'l':
				//Check if another measure was also given as a parameter
				if(msr != none && msr != algnWoutOffs){
					cerr << "ERROR: Only one measure can be calculated at a time!" << endl;

					return false;
				}

				msr = algnWoutOffs;
				break;
			case 'h':
				return false;
			default:
				break;
		}
	}

	return seqaGvn && seqbGvn && msr != none;
}
