#include <parasail.h>
#include <iostream>
#include <cstring>

#include "Alignment.cpp"

int main(){
	//Some fixed alignment parameters
	const int32_t d = 2;
	const int32_t e = 1;
	const parasail_matrix_t *subMat = parasail_matrix_create("ACGT", 1, -1);

	/*IO sequences*/
	const char* query = "ACGTATAGT";
	const int32_t qLen = strlen(query);
	const char* target = "TTTTTTACCTACGTAATTTT";
	const int32_t tLen = strlen(target);
	/*
	typedef struct parasail_string {
	    size_t l;
	    char *s;
	} parasail_string_t;

	typedef struct parasail_sequence {
	    parasail_string_t name;
	    parasail_string_t comment;
	    parasail_string_t seq;
	    parasail_string_t qual;
	} parasail_sequence_t;

	typedef struct parasail_sequences {
	    parasail_sequence_t *seqs;
	    size_t l;
	    size_t characters;
	    size_t shortest;
	    size_t longest;
	    float mean;
	    float stddev;
	} parasail_sequences_t;

	parasail_sequences_t* parasail_sequences_from_file(const char *fname);

	void parasail_sequences_free(parasail_sequences_t *sequences);
	*/

	/* TODO: Make profile from query */

	//Calculate alignment
	parasail_result_t* res = parasail_sg_dx_trace_striped_16(query, qLen, target, tLen, d, e, subMat);

	// int32_t algndBpQ = *(parasail_result_get_length_row(res));
	// int32_t algndBpT = *(parasail_result_get_length_col(res));
	int32_t s = parasail_result_get_score(res);
	int32_t end_q = parasail_result_get_end_query(res);
	int32_t end_t = parasail_result_get_end_ref(res);
	// int32_t alen = parasail_result_get_length(res);
	parasail_cigar_t* cig = parasail_result_get_cigar(res, query, qLen, target, tLen, subMat);
	char* cigStr = parasail_cigar_decode(cig);

	//Testing
	// cout << "sizeof(int): " << sizeof(int) << endl;
	// cout << "sizeof(int32_t): " << sizeof(int32_t) << endl;
	string st = "123";
	cout << "st is " << (st.empty() ? "" : "not ") << "empty" << endl;
	uint32_t i = stoul(st);
	cout << i << endl;

	// typedef struct parasail_cigar_ {
	//     uint32_t *seq;
	//     int len;
	//     int beg_query;
	//     int beg_ref;
	// } parasail_cigar_t;

	cout << "Score: " << s << endl << "End in query: " << end_q << endl << "End in target: " << end_t << endl;
	// cout << "Aligned bases in query: " << *(parasail_result_get_length_row(res)) << endl << "Aligned bases in target: " << *(parasail_result_get_length_row(res)) << endl;
	// cout << "Start in query; " << rCig->beg_query << " Start in target: " << rCig->beg_ref << endl;
	cout << "Cigar string: " << cigStr << endl;

	uint32_t start, end;
	fndMtchngCrds(cigStr, 2 * cig->len, start, end);

	cout << "First match on target at position " << start << " last at " << end << endl;

	start = 0, end = 0;
	uint32_t aLen;
	uint32_t nMtchs;
	prsCgr(*cig, start, end, aLen, nMtchs);

	cout << "start: " << start  << endl;

	parasail_traceback_generic(query, qLen, target, tLen, "Query:", "Target:", subMat, res, '|', '*', 'X', 60, 7, 1);
	//Free alignment results
	parasail_result_free(res);
	res = parasail_sg_dx_stats_striped_16(query, qLen, target, tLen, d, e, subMat);
	cout << "Number of matches is " << parasail_result_get_matches(res) << endl;

	// //Free substitution matrix
	// parasail_matrix_free(subMat);

	return 0;
}