#include <getopt.h>
#include <fstream>

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

//This function parses the program parameters. Returns false if given arguments are not valid
const bool prsArgs(int& nArgs, char** argList, string& pat, string& txt, uint32_t& k, double& hFrac, uint32_t& cw, uint32_t& uw, int32_t& tThres, 
	bool& norm){
	int option_index = 0, a;

	static struct option long_options[] = {
        {"pattern",            required_argument,  0, 'p'},
        {"text",               required_argument,  0, 's'},
        {"ksize",              required_argument,  0, 'k'},
        {"hashratio",          required_argument,  0, 'r'},
        {"commonhashweight",   required_argument,  0, 'c'},
        {"uniquehashweight",   required_argument,  0, 'u'},
        {"hom_thres",          required_argument,  0, 't'},
        {"normalize",          no_argument,        0, 'n'},
        {"help",               no_argument,        0, 'h'},
        {0,                    0,                  0,  0 }
    };

    //Parse all parameters given
	while ((a = getopt_long(nArgs, argList, T_HOM_OPTIONS, long_options, &option_index)) != -1){
		//Assign parameter values
		switch(a){
			case 'p':
				//Save input sequence
				pat = optarg;
				break;
			case 's':
				//Save input sequence
				txt = optarg;
				break;
			case 'k':
				//A k-mer length should be positive
				if(atoi(optarg) <= 0){
					cerr << "ERROR: K-mer length not applicable" << endl;
					return false;
				}

				k = atoi(optarg);
				break;
			case 'r':
				//Check if given value is reasonable to represent a ratio
				if(atof(optarg) <= 0 || atof(optarg) > MAX_RATIO){
					cerr << "ERROR: Given hash ratio not applicable" << endl;

					return false;
				}

				hFrac = atof(optarg);
				break;
			case 'c':
				//Weights should be positive
				if(atoi(optarg) <= 0){
					cerr << "ERROR: Common hash weight not applicable" << endl;

					return false;
				}

				cw = atoi(optarg);
				break;
			case 'u':
				//Weights should be positive
				if(atoi(optarg) <= 0){
					cerr << "ERROR: Unique hash weight not applicable" << endl;

					return false;
				}

				uw = atoi(optarg);
				break;
			case 't':
				tThres = atoi(optarg);
				break;
			case 'n':
				norm = true;
				break;
			case 'h':
				return false;
			default:
				break;
		}
	}

	return !pat.empty() && !txt.empty();
}

//This function reads a file in FASTA format and returns true on success
const bool readFASTA(const string& filePath, string& seq){
	bool headerRead = false;
	string line;

	//Open the file
	ifstream fStr(filePath);

	//Check if the file is open
	if(!fStr.is_open()) return false;

	//Read in queries (one per line)
	while(getline(fStr, line)){
		//We are done if we find a second header in the file
		if(line.front() == '>' && headerRead) break;

		//Header lines are skipped
		if(line.front() == '>'){
			headerRead = true;
			continue;
		}

		
	}
	//Close file
	fStr.close();

	return true;
}
