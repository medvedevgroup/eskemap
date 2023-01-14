#include <iostream>

#include <edlib.h>
#include <parasail/io.h>

using namespace std;

int main(int argc, char **argv){
	parasail_sequence_t s1, s2;
	parasail_sequences_t *inSeqs1, *inSeqs2;
	EdlibAlignResult result;

	if(argc < 3){
		cerr << "Usage: CalcGlobEditDistance <seq1.fasta> <seq2.fasta>" << endl;
		return -1;
	}

	//Load sequences
	inSeqs1 = parasail_sequences_from_file(argv[1]);
	s1 = inSeqs1->seqs[0];
	inSeqs2 = parasail_sequences_from_file(argv[2]);
	s2 = inSeqs2->seqs[0];
	result = edlibAlign(s1.seq.s, s1.seq.l, s2.seq.s, s2.seq.l, edlibDefaultAlignConfig());

	if (result.status == EDLIB_STATUS_OK) cout << result.editDistance << endl;

	edlibFreeAlignResult(result);

	return 0;
}