#include <iostream>
#include <cstring>
#include <vector>
#include <tuple>

// #include <parasail.h>

#include <parasail/io.h>
#include <edlib.h>

// #include "Alignment.cpp"

using namespace std;

int main(int argc, char **argv){
	//Some other fixed constants
	const float THRES = 0.03;

	// const float dec = 0.09552835;
	// const float inter = -73.15283528352302;

	//Input sequences
	uint32_t tPLen;
	char *tPSeq, *tSubP;
	parasail_sequences_t *qrys, *refs;
	parasail_sequence_t q, t;
	//The alignment algorithm's configuration
	EdlibAlignConfig aConf;
	//Piece sequence, sequence length, start in target
	vector<tuple<char*, uint32_t, uint32_t>> targetPieces;
	vector<tuple<char*, uint32_t, uint32_t>>::const_iterator p;
	//Alignment results

	char *cigar;

	uint32_t globOffs, resCounter;
	EdlibAlignResult result;

	/*IO sequences*/
	if(argc != 3){
		cerr << "Usage: FindSimSeqs <Query_FASTA> <Reference_FASTA>" << endl;

		return 1;
	}

	// if(qrys->l > 1) cerr << "WARNING: Multiple query sequences read. Going to process only first one!" << endl;

	//Testing
	// cout << "(int32_t) (THRES * q.seq.l): " << (int32_t) (THRES * q.seq.l) << endl;

	// aConf = edlibNewAlignConfig((int32_t) (THRES * q.seq.l), EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0);
	refs = parasail_sequences_from_file(argv[2]);

	if(refs->l > 1){
		cerr << "ERROR: Multiple reference sequences read. Which one to take?" << endl;

		return 1;
	}

	t = refs->seqs[0];
	qrys = parasail_sequences_from_file(argv[1]);

	//Testing
	// cout << qrys->l << endl;
	// return 0;

	for(size_t i = 0; i < qrys->l; ++i){
		tPSeq = (char*) malloc(t.seq.l + 1);
		strcpy(tPSeq, t.seq.s);
		targetPieces.push_back(make_tuple(tPSeq, t.seq.l, 0));
		q = qrys->seqs[i];
		aConf = edlibNewAlignConfig((int32_t) (THRES * q.seq.l), EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0);
		resCounter = 0;
	
		cout << ">" << q.name.s << " " << q.comment.s << endl;

		// thres = (int32_t) (dec * q.seq.l + inter);
	
		//Testing
		// cout << "thres: " << thres << endl;
		// int32_t a = 2, b = 3;
		// cout << "(float) a / b: " << (float) a / b << endl;
		// return 0;
	
		while(!targetPieces.empty()){
			//Get next target sequence piece
			tPSeq = get<0>(targetPieces.back());
			tPLen = get<1>(targetPieces.back());
			globOffs = get<2>(targetPieces.back());
			targetPieces.pop_back();
			//Calculate alignment
			result = edlibAlign(q.seq.s, q.seq.l, tPSeq, tPLen, aConf);

			//Analyse result
			if(result.status != EDLIB_STATUS_OK){
				cerr << "ERROR: Alignment calculation failed!" << endl;
				return -1;
			}
	
			//Check if result is good enough
			if(result.editDistance > -1){
				//Report result
				cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
	
				cout << globOffs + result.startLocations[0] << " " << globOffs + result.endLocations[0] << " " << cigar << endl;
	
				free(cigar);

				if(++resCounter > 100){
					edlibFreeAlignResult(result);
					free(tPSeq);

					for(p = targetPieces.begin(); p != targetPieces.end(); ++p){
						free(get<0>(*p));
					}

					targetPieces.clear();
					continue;
				}
	
				//Testing
				// cout << "Calculating mapping coordinates..." << endl;
				// cout << globOffs + result.startLocations[0] << " " << globOffs + result.endLocations[0] << endl;
				// cout << "Distance is " << result.editDistance << endl;
	
				//Subdivide target piece
				if(result.startLocations[0] > 0){
					tSubP = (char*) malloc(result.startLocations[0] + 1);
					memcpy(tSubP, &tPSeq[0], result.startLocations[0]);
					tSubP[result.startLocations[0]] = '\0';
					targetPieces.push_back(make_tuple(tSubP, result.startLocations[0], globOffs));
				}
	
				if(result.endLocations[0] < (int32_t) tPLen - 1){
					tSubP = (char*) malloc(tPLen - result.endLocations[0]);
					memcpy(tSubP, &tPSeq[result.endLocations[0] + 1], tPLen - result.endLocations[0] - 1);
					tSubP[tPLen - result.endLocations[0] - 1] = '\0';
					targetPieces.push_back(make_tuple(tSubP, tPLen - result.endLocations[0] - 1, globOffs + result.endLocations[0] + 1));
				}
			}
	
			//Free alignment results
			edlibFreeAlignResult(result);
			free(tPSeq);
		}
	}

	return 0;
}