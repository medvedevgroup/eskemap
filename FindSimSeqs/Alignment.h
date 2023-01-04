#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

//This function parses the given cigar string and calculates the first and last matching character of the first sequence with 
//respect to the second sequence
void fndMtchngCrds(const char* cigStr, const int32_t &strLen, uint32_t &start, uint32_t &end);

//This function parses a parasail cigar and counts the number of leading and trailing deletion, the length of the alignment in bet-
//ween and the total number of matched bases
void prsCgr(const parasail_cigar_t &cgr, uint32_t &start, uint32_t &end, uint32_t &aLen, uint32_t &nMtchs);

#endif