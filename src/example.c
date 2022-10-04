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
	    int kmer_span = iopt.k;

		while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
			mm_idx_stat(mi);

			//Allocate space for first k-mer in indexed sequence
			uint8_t *seq = (uint8_t*) calloc(iopt.k, sizeof(uint8_t));
			//Get sequence of first k-mer
			int seqLen = mm_idx_getseq(mi, 0, 0, iopt.k, seq);

			printf("Length of first %d-mer is %d\nSequence is ", iopt.k, seqLen);
			// for(int i = 0; i < seqLen; ++i){
			// 	printf("%d", seq[i]);
			// }
			// printf("\n");

			//A k-mer we want to query
			// char *str = "TCTTATCTCAGAAAA";
			// char *str = "AAATACACCGCACTG";
			char *str = "TCTTATCTCAGAAAA";
			//2-bit compressed sequence representation of a k-mer
			uint64_t minier = 0;

			for (uint32_t i= 0; i < iopt.k; ++i) {
				int c = seq_nt4_table[(uint8_t)str[i]];
				
				// mm128_t info = { UINT64_MAX, UINT64_MAX };
				
				//Check if current base is ambiguous
				if (c < 4) { // not an ambiguous base
					// int z;
					// if (is_hpc) {
					// 	int skip_len = 1;
					// 	if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					// 		for (skip_len = 2; i + skip_len < len; ++skip_len)
					// 			if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
					// 				break;
					// 		i += skip_len - 1; // put $i at the end of the current homopolymer run
					// 	}
					// 	tq_push(&tq, skip_len);
					// 	kmer_span += skip_len;
					// 	if (tq.count > k) kmer_span -= tq_shift(&tq);
					// } else kmer_span = l + 1 < k? l + 1 : k;

					minier = (minier << 2 | c) & mask;           // forward k-mer
					
					// kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
					
					// if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
					// z = kmer[0] < kmer[1]? 0 : 1; // strand
					// ++l;
					// if (l >= k && kmer_span < 256) {
					// 	info.x = hash64(kmer[z], mask) << 8 | kmer_span;
					// 	info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
					// }
				}// else l = 0, tq.count = tq.front = 0, kmer_span = 0;

				// buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
				// if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
				// 	for (j = buf_pos + 1; j < w; ++j)
				// 		if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
				// 	for (j = 0; j < buf_pos; ++j)
				// 		if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
				// }
				// if (info.x <= min.x) { // a new minimum; then write the old min
				// 	if (l >= w + k && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
				// 	min = info, min_pos = buf_pos;
				// } else if (buf_pos == min_pos) { // old min has moved outside the window
				// 	if (l >= w + k - 1 && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
				// 	for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				// 		if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
				// 	for (j = 0; j <= buf_pos; ++j)
				// 		if (min.x >= buf[j].x) min = buf[j], min_pos = j;
				// 	if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				// 		for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
				// 			if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
				// 		for (j = 0; j <= buf_pos; ++j)
				// 			if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
				// 	}
				// }
				// if (++buf_pos == w) buf_pos = 0;
			}

			printf("minier: %lu\n", minier);
			// minier = 78048;

			//Calculate hash of k-mer to query
	        uint64_t m = hash64(minier, mask);// << 8 | kmer_span;
	        const uint64_t *idx_p = mm_idx_get(mi, m, &n);

	        //Check if k-mer could be found
		    if (n > 0) {
		        //Check if k-mer occurs only once
		        if (n == 1) {
		            //Print sequence id position and strand
		            printf("%ld %d %d", (*idx_p)>>32, ((uint32_t)(*idx_p))>>1, ((uint32_t)*idx_p)&1);
		        }
		        else {
		           	//Iterate over all occurrences
		            for (int i = 0; i < n; i++) {
		                //Print sequence id position and strand
		                printf("%ld %d %d\n", (*idx_p)>>32, ((uint32_t)(*idx_p))>>1, ((uint32_t)*idx_p)&1);
		                //Move to next occurrence
		                idx_p++;
	                }
	            }

	            // printf("m: %lu\n", m);

		            // if(count > 5) break;
		    }

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
		mm_idx_reader_close(r); // close the index reader

		return 0;
	}
