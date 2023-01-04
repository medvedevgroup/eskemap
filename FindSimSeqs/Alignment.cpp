#include <string>

#include "Alignment.h"

using namespace std;

//This function parses the given cigar string and calculates the first and last matching character of the first sequence with 
//respect to the second sequence
void fndMtchngCrds(const char* cigStr, const int32_t &strLen, uint32_t &start, uint32_t &end){
	uint32_t pos = 0;
	string cnt;

	start = 1;
	end = 0;

	//Testing
	cout << "strLen: " << strLen << endl;

	for(int32_t i = 0; i < strLen; ++i){
		switch(cigStr[i]){
			case 'X':
				pos += stoul(cnt);
				cnt.clear();
				break;
			case 'I':
				cnt.clear();
				break;
			case 'D':
				pos += stoul(cnt);
				cnt.clear();
				break;
			case '=':
				if(start > end) start = pos;

				pos += stoul(cnt);
				cnt.clear();

				end = pos - 1;

				//Testing
				cout << end << endl;

				break;
			default:
				cnt.push_back(cigStr[i]);
				break;
		}
	}
}

//This function parses a parasail cigar and counts the number of leading and trailing deletion, the length of the alignment in bet-
//ween and the total number of matched bases
void prsCgr(const parasail_cigar_t &cgr, uint32_t &start, uint32_t &end, uint32_t &aLen, uint32_t &nMtchs){
	//Testing
	// for(uint32_t i = 0; i < cgr.len; ++i){
	// 	cout << (cgr.seq[i] << 28 >> 28) << endl;
	// }

	if(cgr.len > 0 && parasail_cigar_decode_op(cgr.seq[0]) == 'D'){
		start = parasail_cigar_decode_len(cgr.seq[0]);
	}
}
