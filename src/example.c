	// To compile:
	//   gcc -g -O2 example.c libminimap2.a -lz -lm -lpthread

	#include <stdlib.h>
	#include <assert.h>
	#include <stdio.h>
	#include <zlib.h>
	#include "../../software/minimap2/minimap.h"
	#include "../../software/minimap2/kseq.h"
	#include "../../software/minimap2/sketch.c"

	// KSEQ_INIT(gzFile, gzread)

	int main(int argc, char *argv[])
	{
		mm_idxopt_t iopt;
		mm_mapopt_t mopt;
		int n_threads = 1;
		hFrac = 0.2;

		mm_verbose = 2; // disable message output to stderr
		mm_set_opt(0, &iopt, &mopt);
		//Change window size for minimizer calculation
		iopt.w = 1;

		if (argc < 2) {
			fprintf(stderr, "Usage: indexTest <text.fa>\n");
			return 1;
		}

		// open index reader
		mm_idx_reader_t *r = mm_idx_reader_open(argv[1], &iopt, "testindex.idx");
		//The index
		mm_idx_t *mi;
		//Variable to store the number of times a hash is found inside the index
	    int n;
	    //Mask for calculating the hash
	    uint64_t mask = (1ULL<<2*iopt.k) - 1;
	    //I am not super sure what this does...might be relevant for homopolymer compression...
	    int kmer_span = iopt.k, shift1 = 2 * (iopt.k - 1), z;
	    //A variable to store a k-mer's hash
	    uint64_t m;

		while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
			mm_idx_stat(mi);

		    //Allocate space for first k-mer in indexed sequence
		    uint8_t *seq = (uint8_t*) calloc(mi->seq[0].len, sizeof(uint8_t));
			//Get sequence of first k-mer
			int seqLen = mm_idx_getseq(mi, 0, 0, mi->seq[0].len, seq);

			// printf("Length of indexed sequence is %d\n", seqLen);

			//2-bit compressed sequence representation of a k-mer
			uint64_t minier = 0, refkm = 0, revkm = 0;

			for (uint32_t i= 0; i < seqLen; ++i) {
				int c = seq_nt4_table[(uint8_t)seq[i]];
				
				//Check if current base is ambiguous
				if (c < 4) { // not an ambiguous base
					refkm = (refkm << 2 | c) & mask;           // forward k-mer
					revkm = (revkm >> 2) | (3ULL^c) << shift1; // reverse k-mer
				}

				//Check if we have read in the next complete k-mer
				if(i > iopt.k - 2){
					printf("K-mer is %lu\n", refkm);

					//Check if k-mer and its reverse complement are identical
					if(refkm == revkm){
						printf("But we won't find it because it is its own reverse complement");
						continue;
					}

					//Check if reverse complement is smaller
					// if(refkm < revkm){

					//We always want to use the k-mer coming from the reference strand
					minier = refkm;

					// } else{
					// 	minier = revkm;

					// 	printf("But we query its reverse complement\n");
					// }

					//Calculate hash of k-mer to query
					m = hash64(minier, mask);
					//Query hashed k-mer in index
					const uint64_t *idx_p = mm_idx_get(mi, m, &n);

					// printf("Querying done\n");

			        //Check if k-mer could be found
				    if (n > 0) {
				        //Check if k-mer occurs only once
				        if (n == 1) {
				        	// printf("Single occurrence\n");

				            //Print sequence id position and strand
				            printf("%ld %d %d\n", (*idx_p)>>32, ((uint32_t)(*idx_p))>>1, ((uint32_t)*idx_p)&1);
				        }
				        else {
				        	// printf("Multiple occurrences\n");

				           	//Iterate over all occurrences
				            for (int i = 0; i < n; i++) {
				                //Print sequence id position and strand
				                printf("%ld %d %d\n", (*idx_p)>>32, ((uint32_t)(*idx_p))>>1, ((uint32_t)*idx_p)&1);
				                //Move to next occurrence
				                idx_p++;
			                }
			            }
				    }

				    // n = 0;
				}
			}

			// minier = 78048;

			// //Calculate hash of k-mer to query
	  //       uint64_t m = hash64(minier, mask);// << 8 | kmer_span;
	        
	        // printf("hash(TTATCCCAGATCAAC): %lu\n", hash64(minier, mask));
	        // uint8_t count = 0;

	        // for(uint64_t m = 946913546; m < UINT64_MAX; ++m){
	        // 	//Query the index
		       //  const uint64_t *idx_p = mm_idx_get(mi, m, &n);

		       //  //Check if k-mer could be found
		       //  if (n > 1) {
		       //  	//Check if k-mer occurs only once
		       //      if (n == 1) {
		       //      	//Print sequence id position and strand
		       //          printf("%ld %d %d", (*idx_p)>>32, ((uint32_t)(*idx_p))>>1, ((uint32_t)*idx_p)&1);
		       //      }
		       //      else {
		       //      	//Iterate over all occurrences
		       //          for (int i = 0; i < n; i++) {
		       //          	//Print sequence id position and strand
		       //              printf("%ld %d %d\n", (*idx_p)>>32, ((uint32_t)(*idx_p))>>1, ((uint32_t)*idx_p)&1);
		       //              //Move to next occurrence
		       //              idx_p++;
		       //          }
		       //      }

		       //      printf("m: %lu\n", m);

		       //      // if(count > 5) break;
		       //      break;
		       //  }
	        // }

			mm_idx_destroy(mi);
		}
		// mm_idx_reader_close(r); // close the index reader

		return 0;
	}
