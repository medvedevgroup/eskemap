// To compile:
//   gcc -g -O2 example.c libminimap2.a -lz -lm -lpthread

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "../../software/minimap2/minimap.h"
#include "../../software/minimap2/kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	int n_threads = 1;

	mm_verbose = 2; // disable message output to stderr
	mm_set_opt(0, &iopt, &mopt);
	//Change window size for minimizer calculation
	iopt.w = 1;

	// mopt.flag |= MM_F_CIGAR; // perform alignment

	if (argc < 2) {
		fprintf(stderr, "Usage: indexTest <text.fa>\n");
		return 1;
	}

	// // open query file for reading; you may use your favorite FASTA/Q parser
	// gzFile f = gzopen(argv[2], "r");
	// assert(f);
	// kseq_t *ks = kseq_init(f);

	// open index reader
	mm_idx_reader_t *r = mm_idx_reader_open(argv[1], &iopt, "testindex.idx");
	mm_idx_t *mi;
	while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
	// 	mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
	// 	mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
	// 	gzrewind(f);
	// 	kseq_rewind(ks);
	// 	while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
	// 		mm_reg1_t *reg;
	// 		int j, i, n_reg;
	// 		reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query
	// 		for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
	// 			mm_reg1_t *r = &reg[j];
	// 			assert(r->p); // with MM_F_CIGAR, this should not be NULL
	// 			printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
	// 			printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:", mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq);
	// 			for (i = 0; i < r->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
	// 				printf("%d%c", r->p->cigar[i]>>4, MM_CIGAR_STR[r->p->cigar[i]&0xf]);
	// 			putchar('\n');
	// 			free(r->p);
	// 		}
	// 		free(reg);
		// }
	// 	mm_tbuf_destroy(tbuf);

		// int mm_idx_name2id(const mm_idx_t *mi, const char *name);
		// printf("w is %d\n", iopt.w);
		mm_idx_stat(mi);
		uint8_t *seq = (uint8_t*) calloc(mi->seq[0].len + 1, sizeof(uint8_t));
		int seqLen = mm_idx_getseq(mi, 0, 0, mi->seq[0].len, seq);

		for(int i = 0; i < seqLen; ++i){
			printf("%d ", seq[i]);
		}
		// printf("seqLen: %d\n", seqLen);
		
		// int mm_idx_alt_read(mm_idx_t *mi, const char *fn);
		// int mm_idx_bed_read(mm_idx_t *mi, const char *fn, int read_junc);
		// int mm_idx_bed_junc(const mm_idx_t *mi, int32_t ctg, int32_t st, int32_t en, uint8_t *s);
		mm_idx_destroy(mi);
	}
	mm_idx_reader_close(r); // close the index reader
	// kseq_destroy(ks); // close the query file
	// gzclose(f);
	return 0;
}
