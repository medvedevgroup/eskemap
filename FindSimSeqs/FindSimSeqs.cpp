#include <iostream>
#include <cstring>
#include <vector>
#include <tuple>

#include <parasail.h>
#include <parasail/io.h>

#include "Alignment.cpp"

int main(int argc, char **argv){
	//Some fixed alignment parameters
	const int32_t d = 2;
	const int32_t e = 1;
	const parasail_matrix_t *subMat = parasail_matrix_create("ACGT", 1, -1);
	//Some other fixed constants
	const float MIN_SID = 0.95;
	//Input sequences
	uint32_t tPLen;
	char *tPSeq, *tSubP;
	parasail_sequences_t *qrys, *refs;
	parasail_sequence_t q, t;
	//Piece sequence, sequence length, start in target
	vector<tuple<char*, uint32_t, uint32_t>> targetPieces;
	//Alignment results
	uint32_t start, end, globOffs;
	float aLen, nMtchs;
	parasail_result_t* res;
	parasail_cigar_t* cig;

	/*IO sequences*/

	//Testing
	// const char* query = "ACGTATAGT";
	// const char* query = "CACATAACGACGCCTCGCATCGCGTTACGCTCACGGCCAGAGGCTGCAG";
	// const int32_t qLen = strlen(query);
	// // const char* target = "TTTTTTACCTACGTAATTTT";
	// const char* target = "TTTCTTGCCGATCACATCCGACGCGTCTCGCAAACGGGTTACGCTCAAGGCAGCATAACGGCCAGAGGCTGCAG";
	// const int32_t tLen = strlen(target);

	if(argc != 3){
		cerr << "Usage: FindSimSeqs <Query_FASTA> <Reference_FASTA>" << endl;

		return 1;
	}

	//Testing
	// cout << "argv[1]: " << argv[1] << endl;

	qrys = parasail_sequences_from_file(argv[1]);

	//Testing
	// qrys = parasail_sequences_from_file("../../simulations/randSeqCopy_l50000_rid0_m0.005_d0.01_i0.01_cn2.fasta");

	if(qrys->l > 1) cerr << "WARNING: Multiple query sequences read. Going to process only first one!" << endl;

	q = qrys->seqs[0];
	refs = parasail_sequences_from_file(argv[2]);

	//Testing
	// refs = parasail_sequences_from_file("../../simulations/randSeqCopy_l50000_rid0_m0.005_d0.01_i0.01_cn2.fasta");

	if(refs->l > 1){
		cerr << "ERROR: Multiple reference sequences read. Which one to take?" << endl;

		return 1;
	}

	t = refs->seqs[0];
	targetPieces.push_back(make_tuple(t.seq.s, t.seq.l, 0));

	//This causes trouble
	// parasail_sequences_free(refs);

	//Testing
	// cout << "Query sequence: " << q.seq.s << endl;
	// cout << "Query length: " << q.seq.l << endl;
	// // return 0;
	// cout << "Target sequence: " << t.seq.s << endl;
	// // return 0;
	// cout << "Target length: " << t.seq.l << endl;

	/* TODO: Make profile from query */
	//no idea how this should work...
	// qP = parasail_profile_create_16 (q.seq.s, q.seq.l, subMat);

	while(!targetPieces.empty()){
		//Get next target sequence piece
		tPSeq = get<0>(targetPieces.back());
		tPLen = get<1>(targetPieces.back());
		globOffs = get<2>(targetPieces.back());
		targetPieces.pop_back();
		//Calculate alignment
		res = parasail_sg_dx_trace_striped_16(q.seq.s, q.seq.l, tPSeq, tPLen, d, e, subMat);
		//Analyse result
		cig = parasail_result_get_cigar(res, q.seq.s, q.seq.l, tPSeq, tPLen, subMat);
		prsCgr(*cig, start, end, aLen, nMtchs);

		//Testing
		// char* cigStr = parasail_cigar_decode(cig);
		// cout << "Cigar string: " << cigStr << endl;
		// // cout << "targetPieces.size(): " << targetPieces.size() << endl;
		// if(globOffs == 2997){
		// 	cout << "tPSeq: " << tPSeq << endl;
		// 	cout << "tPLen: " << tPLen << endl;
		// }
		// cout << "aLen: " << aLen << endl;
		// cout << "nMtchs: " << nMtchs << endl;

		//Check if result is good enough
		if(nMtchs / max(aLen, (float) q.seq.l) >= MIN_SID){
			//Report result
			cout << globOffs + start << " " << globOffs + end << " " << parasail_cigar_decode(cig) << endl;

			//Subdivide target piece

			//Testing
			// cout << "tPLen: " << tPLen << endl;
			// cout << "globOffs: " << globOffs << endl;
			// cout << "end: " << end << endl;
			// cout << "start: " << start << endl;

			if(start > 0){
				tSubP = (char*) malloc(start + 1);
				memcpy(tSubP, &tPSeq[0], start);
				tSubP[start] = '\0';
				targetPieces.push_back(make_tuple(tSubP, start, globOffs));

				//Testing
				// cout << "tSubP: " << tSubP << endl;
				// char *a = "T";
				// targetPieces.push_back(make_tuple(a, 1, 0));
			}

			if(end < tPLen - 1){
				tSubP = (char*) malloc(tPLen - end);
				memcpy(tSubP, &tPSeq[end + 1], tPLen - end - 1);
				tSubP[tPLen - end - 1] = '\0';
				targetPieces.push_back(make_tuple(tSubP, tPLen - end - 1, globOffs + end + 1));

				//Testing
				// char *b = "A";
				// targetPieces.push_back(make_tuple(b, 1, 21));
				// cout << "tSubP: " << tSubP << endl;
			}
		}

		//Testing
		// cout << "Sequence identity: " << nMtchs / max(aLen, (float) q.seq.l) << endl;
		// if(globOffs == 2997){
		// 	cout << "globOffs: " << globOffs << endl;
		// parasail_traceback_generic(q.seq.s, q.seq.l, tPSeq, tPLen, "Query:", "Target:", subMat, res, '|', '*', 'X', 60, 7, 1);
		// }

		// free(tPSeq);
		parasail_cigar_free(cig);
		//Free alignment results
		parasail_result_free(res);
	}

	//Testing
	// return 0;
	// int32_t algndBpQ = *(parasail_result_get_length_row(res));
	// int32_t algndBpT = *(parasail_result_get_length_col(res));
	// int32_t s = parasail_result_get_score(res);
	// int32_t end_q = parasail_result_get_end_query(res);
	// int32_t end_t = parasail_result_get_end_ref(res);
	// int32_t alen = parasail_result_get_length(res);
	// cout << "sizeof(int): " << sizeof(int) << endl;
	// cout << "sizeof(int32_t): " << sizeof(int32_t) << endl;
	// string st = "123";
	// cout << "st is " << (st.empty() ? "" : "not ") << "empty" << endl;
	// uint32_t i = stoul(st);
	// cout << i << endl;
	// cout << "Score: " << s << endl << "End in query: " << end_q << endl << "End in target: " << end_t << endl;
	// cout << "Aligned bases in query: " << *(parasail_result_get_length_row(res)) << endl << "Aligned bases in target: " << *(parasail_result_get_length_row(res)) << endl;
	// cout << "Start in query; " << rCig->beg_query << " Start in target: " << rCig->beg_ref << endl;
	// fndMtchngCrds(cigStr, 2 * cig->len, start, end);
	// cout << "First match on target at position " << start << " last at " << end << endl;
	// cout << "Start: " << start << " end: " << end << " aLen: " << aLen << " nMtchs: " << nMtchs << endl;
	// res = parasail_sg_dx_stats_striped_16(q.seq.s, q.seq.l, t.seq.s, t.seq.l, d, e, subMat);
	// cout << "Number of matches is " << parasail_result_get_matches(res) << endl;

	// //Free substitution matrix
	// parasail_matrix_free(subMat);

	return 0;
}