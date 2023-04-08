#include <iostream>

#include <edlib.h>
#include <parasail/io.h>

using namespace std;

int main(int argc, char **argv){
	// parasail_sequence_t s1, s2;
	parasail_sequences_t *inSeqs1, *inSeqs2;
	EdlibAlignResult result;

	if(argc < 3){
		cerr << "Usage: CalcGlobEditDistance <seq1.fasta> <seq2.fasta>" << endl;
		return -1;
	}

	//Load sequences
	inSeqs1 = parasail_sequences_from_file(argv[1]);
	inSeqs2 = parasail_sequences_from_file(argv[2]);

	//Iterate over all pairs of sequences
	for(uint32_t i = 0; i < inSeqs1->l; ++i){
		parasail_sequence_t s1 = inSeqs1->seqs[i];

		//Testing
		// cout << "Sequence 1: ";
		// for(uint32_t k = 0; k < 10; ++k) cout << inSeqs1->seqs[i].seq.s[k];
		// cout << endl;

		for(uint32_t j = 0; j < inSeqs2->l; ++j){
			parasail_sequence_t s2 = inSeqs2->seqs[j];

			//Testing
			// cout << "Sequence 2: ";
			// for(uint32_t k = 0; k < 10; ++k) cout << inSeqs2->seqs[j].seq.s[k];
			// cout << endl;

			result = edlibAlign(s1.seq.s, s1.seq.l, s2.seq.s, s2.seq.l, edlibDefaultAlignConfig());
			
			if (result.status == EDLIB_STATUS_OK) cout << result.editDistance << endl;
		
			edlibFreeAlignResult(result);
		}
	}

	return 0;
}