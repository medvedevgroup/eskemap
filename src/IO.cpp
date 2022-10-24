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
const bool prsArgs(int& nArgs, char** argList, string& pFl, string& tFl, uint32_t& k, double& hFrac, uint32_t& cw, uint32_t& uw, 
	int32_t& tThres, bool& norm){
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
				pFl = optarg;
				break;
			case 's':
				//Save input sequence
				tFl = optarg;
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

	return !pFl.empty() && !tFl.empty();
}

//This function reads a file in FASTA format and returns true on success
//NOTE: Using this function, we first read in a complete sequence before calculating a sketch from it. Calculating a sketch could
//		also be done on the fly while reading the file. This possibility would be more memory saving most likely, because we would
//		not have to store the whole sequence in memory first. However, we had no chance to estimate the sketch size. Thus, the vec-
//		tor to store the sketch would have to be elongated several times which might be time consuming. On the other hand, it might
//		also be time consuming to handle the complete sequence in memory first. Eventually, it will depend on a try to find out
//		what is the better alternative. Could also be that it does not matter at all!
const bool readFASTA(const string& filePath, string& seq){
	bool headerRead = false, lnBrkDiscvd = false;
	char c;

	//Open the file
	ifstream fStr(filePath);

	//Check if the file is open
	if(!fStr.is_open()) return false;

	//Read in file character by character
	while(fStr.get(c)){
		//We are done if we find a second header (which can only start after at least one line break) in the file
		if(c == '>' && headerRead && lnBrkDiscvd) break;

		//Header lines are skipped
		if(c == '>'){
			headerRead = true;
			continue;
		}

		//The first line break indicates that we have left the header line
		if(c == '\n'){
			lnBrkDiscvd = true;
			continue;
		}

		//There is no sequence to load in the header line
		if(headerRead && !lnBrkDiscvd) continue;

		//We are only interested in unambigous, unmasked nucleotides
		if(c == 'A' || c == 'C' || c == 'G' || c == 'T') seq += c;
	}

	//Close file
	fStr.close();

	return true;
}

//This function reads in batches of FASTA sequence entries from file and transforms them into sketches. Returns false if end of file
//was reached.
const bool lPttnSks(ifstream& fStr, const uint32_t& k, const double& hFrac, vector<pair<string, Sketch>>& pSks){
	bool headerRead, idRead = false, lnBrkDiscvd = false;
	char c;
	string seq, seqID;

	//Check if the file is open
	if(!fStr.is_open()){
		//Testing
		// cout << "1 Option 2" << endl;

		return false;
	}

	//Testing
	// cout << "1 Option 1" << endl;

	//If this function is called iteratively, the '>' of the next entry has already been read
	headerRead = fStr.gcount() != 0;

	//Read in file character by character
	while(fStr.get(c)){
		//Testing
		// cout << "lPttnSks: Start of while loop" << endl;

		//An entry's sequence is complete if we find a second header (which can only start after at least one line break) in the 
		//file
		if(c == '>' && headerRead && lnBrkDiscvd){
			//Testing
			// cout << "2 Option 1" << endl;
			// exit(0);

			//Add sequence's sketch and id to result vector
			pSks.push_back(make_pair(seqID, buildSketch(seq, k, hFrac)));
			//Clear sequence id
			seqID.clear();
			//Clear sequence
			seq.clear();

			//Check if enough sequences have been read
			if(pSks.size() == PATTERN_BATCH_SIZE){
				//Testing
				// cout << "3 Option 1" << endl;

				return true;
			}

			//Testing
			// cout << "3 Option 2" << endl;

			//Reset id-read flag
			idRead = false;
			//Reset line-break-discovered flag
			lnBrkDiscvd = false;
			continue;
		}

		//Testing
		// cout << "2 Option 2" << endl;
		// cout << "5 Option " << (headerRead ? "1" : "2") << endl;
		// cout << "6 Option " << (!headerRead && c == '>' ? "1" : "2") << endl;
		// cout << "7 Option " << (idRead ? "1" : "2") << endl;
		// cout << "8 Option " << (!idRead && (headerRead && c == ' ' && !lnBrkDiscvd) ? "1" : "2") << endl;
		// cout << "9 Option " << (lnBrkDiscvd ? "1" : "2") << endl;
		// cout << "10 Option " << (!lnBrkDiscvd && c == '\n' ? "1" : "2") << endl;
		
		//Note if we have completely read the sequence id
		idRead = idRead || (headerRead && c == ' ' && !lnBrkDiscvd);
		//Note if we have found the first line break after a new header started
		lnBrkDiscvd = lnBrkDiscvd || c == '\n';

		//Update sequence id if we are still reading it
		if(headerRead && !lnBrkDiscvd && !idRead){
			//Testing
			// cout << "4 Option 1" << endl;

			seqID += c;
			continue;
		}

		//Note if we have found the beginning of a header
		headerRead = headerRead || (c == '>');

		//Testing
		// cout << "4 Option 2" << endl;

		//There is no sequence to load in the header line
		if(headerRead && !lnBrkDiscvd){
			//Testing
			// cout << "11 Option 1" << endl;

			continue;
		}

		//Testing
		// cout << "11 Option 2" << endl;

		//We are only interested in unambigous, unmasked nucleotides
		if(c == 'A' || c == 'C' || c == 'G' || c == 'T'){
			//Testing
			// cout << "12 Option 1" << endl;

			seq += c;
		} else{
			//Testing
			// cout << "12 Option 2" << endl;
		}
	}

	//Testing
	// cout << "lPttnSks: seq: " << seq << endl;
	// exit(0);

	//Add last entry's sketch and sequence id to result vector if it is not empty
	if(!seq.empty()){
		//Testing
		// cout << "13 Option 1" << endl;

		pSks.push_back(make_pair(seqID, buildSketch(seq, k, hFrac)));
	} else{
		//Testing
		// cout << "13 Option 2" << endl;
	}

	return false;
}
