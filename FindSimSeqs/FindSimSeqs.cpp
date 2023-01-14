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
	// Some fixed alignment parars
	// const int32_t d = 2;
	// const int32_t e = 1;
	// const parasail_matrix_t *subMat = parasail_matrix_create("ACGT", 1, -1);

	//Some other fixed constants

	// const float MIN_SID = 0.95;

	const float dec = 0.09552835;
	const float inter = -73.15283528352302;
	//Input sequences
	uint32_t tPLen;
	char *tPSeq, *tSubP;//, *tPsubSeq;
	parasail_sequences_t *qrys, *refs;
	parasail_sequence_t q, t;
	//Piece sequence, sequence length, start in target
	vector<tuple<char*, uint32_t, uint32_t>> targetPieces;
	//Alignment results
	char *cigar;
	uint32_t globOffs;//start, end, subSeqLen, , aEnd
	int32_t thres;
	EdlibAlignResult result;

	// float aLen, nMtchs;
	// parasail_result_t* res;
	// parasail_cigar_t* cig;

	/*IO sequences*/
	if(argc != 3){
		cerr << "Usage: FindSimSeqs <Query_FASTA> <Reference_FASTA>" << endl;

		return 1;
	}

	qrys = parasail_sequences_from_file(argv[1]);

	if(qrys->l > 1) cerr << "WARNING: Multiple query sequences read. Going to process only first one!" << endl;

	q = qrys->seqs[0];
	refs = parasail_sequences_from_file(argv[2]);

	if(refs->l > 1){
		cerr << "ERROR: Multiple reference sequences read. Which one to take?" << endl;

		return 1;
	}

	t = refs->seqs[0];
	targetPieces.push_back(make_tuple(t.seq.s, t.seq.l, 0));

	//Testing
	// result = edlibAlign(q.seq.s, q.seq.l, t.seq.s, t.seq.l, edlibNewAlignConfig(1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
	// if (result.status == EDLIB_STATUS_OK){
	// 	cout << "Edit distance: " << result.editDistance << endl;
	// 	for(int32_t i = 0; i < result.numLocations; ++i){
	// 		cout << "Start in Reference: " << result.startLocations[i] << endl;
	// 		cout << "End in Reference: " << result.endLocations[i] << endl;
	// 	}
	// 	cout << "Alignment length: " << result.alignmentLength << endl;
	// 	cout << "Query: " << q.seq.s << endl;
	// 	cout << "Reference: " << t.seq.s << endl;
	// 	cout << "Cigar: " << edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD) << endl;
	// }
	// edlibFreeAlignResult(result);
	// return 0;
	// int32_t a = 32;
	// cout << "a: " << a << endl;
	// cout << "(int32_t) (-0.13 * a - 17.84): " << (int32_t) (-0.13 * a - 17.84) << endl;
	// return 0;

	thres = (int32_t) (dec * q.seq.l + inter);

	//Testing
	cout << "thres: " << thres << endl;

	while(!targetPieces.empty()){
		//Get next target sequence piece
		tPSeq = get<0>(targetPieces.back());
		tPLen = get<1>(targetPieces.back());
		globOffs = get<2>(targetPieces.back());
		targetPieces.pop_back();
		//Calculate alignment
		result = edlibAlign(q.seq.s, q.seq.l, tPSeq, tPLen, edlibNewAlignConfig(thres, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

		//Analyse result
		if(result.status != EDLIB_STATUS_OK){
			cerr << "ERROR: Alignment calculation failed!" << endl;
			return -1;
		}

		// aEnd = parasail_result_get_end_ref(res);
		// subSeqLen = min(2 * parasail_result_get_length(res), (int32_t) aEnd + 1);
		// tPsubSeq = (char*) malloc(subSeqLen + 1);

		// //Testing
		// // cout << "(aEnd + 1) - subSeqLen: " << (aEnd + 1) - subSeqLen << endl;
		// // cout << "tPSeq[(aEnd + 1) - subSeqLen]: " << tPSeq[(aEnd + 1) - subSeqLen] << endl;

		// memcpy(tPsubSeq, &tPSeq[(aEnd + 1) - subSeqLen], subSeqLen);
		// tPsubSeq[subSeqLen] = '\0';

		// //Testing
		// // cout << "tPsubSeq: " << tPsubSeq << endl;
		
		// parasail_result_free(res);
		// res = parasail_sg_dx_trace_striped_sat(q.seq.s, q.seq.l, tPsubSeq, subSeqLen, d, e, subMat);
		// cig = parasail_result_get_cigar(res, q.seq.s, q.seq.l, tPsubSeq, subSeqLen, subMat);
		// prsCgr(*cig, start, end, aLen, nMtchs);

		// if(subSeqLen < aEnd + 1){
		// 	start += (aEnd + 1) - subSeqLen;
		// 	end += (aEnd + 1) - subSeqLen;
		// }

		//Check if result is good enough
		if(result.editDistance <= thres){
			//Report result
			cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);

			cout << globOffs + result.startLocations[0] << " " << globOffs + result.endLocations[0] << " " << cigar << endl;

			free(cigar);

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

		// free(tPsubSeq);

		free(tPSeq);
	}

	return 0;
}