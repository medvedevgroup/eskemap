#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <unordered_map>
#define __STDC_LIMIT_MACROS
#include "kvec.h"
#include "mmpriv.h"

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

typedef struct { // a simplified version of kdq
	int front, count;
	int a[32];
} tiny_queue_t;

static inline void tq_push(tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & 0x1f] = x;
}

static inline int tq_shift(tiny_queue_t *q)
{
	int x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= 0x1f;
	--q->count;
	return x;
}

// double hFrac = 0.1;

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
// void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
// {
// 	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0}, hcnt = 0;
// 	//Calculate maximum hash value in sketch
// 	const uint64_t maxHash = hFrac * mask;
// 	int i, l, buf_pos, min_pos, kmer_span = 0;//, j
// 	mm128_t buf[256];//, min = { UINT64_MAX, UINT64_MAX };
// 	tiny_queue_t tq;

// 	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
// 	memset(buf, 0xff, w * 16);
// 	memset(&tq, 0, sizeof(tiny_queue_t));
// 	kv_resize(mm128_t, km, *p, p->n + len/w);

// 	//Testing
// 	// printf("mm_sketch: Sequence to index: %s\n", str);
// 	// printf("mm_sketch: Iterating over text sketch\n");

// 	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
// 		int c = seq_nt4_table[(uint8_t)str[i]];
// 		mm128_t info = { UINT64_MAX, UINT64_MAX };
// 		//We always consider the reference strand
// 		int z = 0;

// 		if (c < 4) { // not an ambiguous base
// 			// if (is_hpc) {
// 			// 	int skip_len = 1;
// 			// 	if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
// 			// 		for (skip_len = 2; i + skip_len < len; ++skip_len)
// 			// 			if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
// 			// 				break;
// 			// 		i += skip_len - 1; // put $i at the end of the current homopolymer run
// 			// 	}
// 			// 	tq_push(&tq, skip_len);
// 			// 	kmer_span += skip_len;
// 			// 	if (tq.count > k) kmer_span -= tq_shift(&tq);
// 			// } else 
// 			kmer_span = l + 1 < k? l + 1 : k;
// 			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
// 			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer

// 			//We do not skip anything...
// 			// if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand

// 			//Testing
// 			// printf("K-mer hash: %lu hcnt: %lu\n", hash64(kmer[z], mask), hcnt);

// 			// z = kmer[0] < kmer[1]? 0 : 1; // strand
// 			++l;
// 			if (l >= k && kmer_span < 256) {
// 				//Well, we want to skip blacklisted k-mers
// 				if(bLstmers.contains(hash64(kmer[z], mask))) continue;

// 				info.x = hash64(kmer[z], mask) << 8 | kmer_span;
// 				info.y = (uint64_t)rid<<32 | (uint32_t)hcnt<<1 | z;
// 			}
// 		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;

// 		//Testing
// 		// if(hash64(kmer[z], mask) == 49938615) printf("mm_sketch: Found hash %lu in text sketch\n", hash64(kmer[z], mask));

// 		//Check if k-mers hash is small enough
// 		if(l >= k && hash64(kmer[z], mask) <= maxHash){
// 			//Testing
// 			// if(hash64(kmer[z], mask) == 49938615) printf("mm_sketch: The hash is pushed to the index\n");
// 			// printf("%lu ", hash64(kmer[z], mask));

// 			//Store k-mer inside index
// 			kv_push(mm128_t, km, *p, info);
// 			++hcnt;
// 		}

// 		// buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
// 		// if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
// 		// 	for (j = buf_pos + 1; j < w; ++j)
// 		// 		if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
// 		// 	for (j = 0; j < buf_pos; ++j)
// 		// 		if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
// 		// }
// 		// if (info.x <= min.x) { // a new minimum; then write the old min
// 		// 	if (l >= w + k && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
// 		// 	min = info, min_pos = buf_pos;
// 		// } else if (buf_pos == min_pos) { // old min has moved outside the window
// 		// 	if (l >= w + k - 1 && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
// 		// 	for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
// 		// 		if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
// 		// 	for (j = 0; j <= buf_pos; ++j)
// 		// 		if (min.x >= buf[j].x) min = buf[j], min_pos = j;
// 		// 	if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
// 		// 		for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
// 		// 			if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
// 		// 		for (j = 0; j <= buf_pos; ++j)
// 		// 			if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
// 		// 	}
// 		// }
// 		// if (++buf_pos == w) buf_pos = 0;
// 	}

// 	//Testing
// 	// printf("\n");
// 	// printf("mm_sketch: maxHash: %lu\n", maxHash);

// 	// if (min.x != UINT64_MAX)
// 	// 	kv_push(mm128_t, km, *p, min);
// }

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0}, hcnt = 0;
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	memset(buf, 0xff, w * 16);
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	//Testing
	// printf("mm_sketch: len: %d\n", len);
	// printf("mm_sketch: sizeof(int): %lu\n", sizeof(int));

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info = { UINT64_MAX, UINT64_MAX };
		if (c < 4) { // not an ambiguous base
			int z;
			if (is_hpc) {
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k) kmer_span -= tq_shift(&tq);
			} else 
			kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				// Well, we want to skip blacklisted k-mers
				// if(bLstmers.contains(hash64(kmer[z], mask))) continue;

				//Testing                  232168040628
				// if(hash64(kmer[z], mask) == 1389627116){
				// // 	printf("Hash was calculated\n");
				// // 	printf("Hash for reference strand was %lu\n", hash64(kmer[0], mask));
				// // 	printf("z: %d\n", z);
				// 	printf("i:%d\n", i);
				// }
				// if(hash64(kmer[0], mask) == 47147373989) printf("The corresponding hash on the reference strand was here\n");

				info.x = hash64(kmer[z], mask) << 8 | kmer_span;
				// info.y = (uint64_t)rid<<32 | (uint32_t)hcnt<<1 | z;
				info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;
			}
		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j){
				if (min.x == buf[j].x && buf[j].y != min.y && !bLstmers.contains(buf[j].x >> 8)){
					//Testing
					// if((buf[j].x >> 8) == 1389627116){
					// 	printf("i:%d\n", i);
					// }
					// if(i >= 62396600 && i <=62396619){
					// printf("Hash in iteration i=%d: %lu\n", i, buf[j].x >> 8);
					// }
					//Print hash
					// printf("%lu\n", buf[j].x >> 8);
					//Print position
					// printf("%lu\n", ((buf[j].y << 32) >> 33) + 1 - k);

					int zp = (int) ((buf[j].y << 63) >> 63);
					buf[j].y = (buf[j].y >> 32) << 32 | (uint32_t)hcnt<<1 | zp;
					kv_push(mm128_t, km, *p, buf[j]);
					++hcnt;
				}
			}
			for (j = 0; j < buf_pos; ++j){
				if (min.x == buf[j].x && buf[j].y != min.y && !bLstmers.contains(buf[j].x >> 8)){
					//Testing
					// if((buf[j].x >> 8) == 1389627116){
					// 	printf("i:%d\n", i);
					// // printf("%lu\n", buf[j].x >> 8);
					// // 	printf("second push\n");
					// }
					// if(i >= 62396600 && i <=62396619){
					// printf("Hash in iteration i=%d: %lu\n", i, buf[j].x >> 8);
					// }
					//Print hash
					// printf("%lu\n", buf[j].x >> 8);
					//Print position
					// printf("%lu\n", ((buf[j].y << 32) >> 33) + 1 - k);

					int zp = (int) ((buf[j].y << 63) >> 63);
					buf[j].y = (buf[j].y >> 32) << 32 | (uint32_t)hcnt<<1 | zp;
					kv_push(mm128_t, km, *p, buf[j]);
					++hcnt;
				}
			}
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX && !bLstmers.contains(min.x >> 8)){
				//Testing          2567431378
				// if((min.x >> 8) == 1389627116){
				// 	printf("i:%d\n", i);
				// // 	printf("Third push\n");
				// // 	printf("info.x: %lu\n", info.x >> 8);
				// }
				// if(i >= 62396600 && i <=62396619){
				// printf("Hash in iteration i=%d: %lu\n", i, min.x >> 8);
				// }
				//Print hash
				// printf("%lu\n", min.x >> 8);
				//Print position
				// printf("%lu\n", ((min.y << 32) >> 33) + 1 - k);

				int zp = (int) ((min.y << 63) >> 63);
				min.y = (min.y >> 32) << 32 | (uint32_t)hcnt<<1 | zp;
				kv_push(mm128_t, km, *p, min);
				++hcnt;
			}
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1 && min.x != UINT64_MAX && !bLstmers.contains(min.x >> 8)){
				//Testing
				// if((min.x >> 8) == 1389627116){
				// 	printf("i:%d\n", i);
				// 	// printf("Fourth push\n");
				// }
				// if(i >= 62396600 && i <=62396619){
				// printf("Hash in iteration i=%d: %lu\n", i, min.x >> 8);
				// }
				//Print hash
				// printf("%lu\n", min.x >> 8);
				//Print position
				// printf("%lu\n", ((min.y << 32) >> 33) + 1 - k);

				int zp = (int) ((min.y << 63) >> 63);
				min.y = (min.y >> 32) << 32 | (uint32_t)hcnt<<1 | zp;
				kv_push(mm128_t, km, *p, min);
				++hcnt;
			}
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j){ // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y && !bLstmers.contains(buf[j].x >> 8)){
						//Testing
						// if((buf[j].x >> 8) == 1389627116){
						// 	printf("i:%d\n", i);
						// // 	printf("Fiveth push\n");
						// }
						// if(i >= 62396600 && i <=62396619){
						// printf("Hash in iteration i=%d: %lu\n", i, buf[j].x >> 8);
						// }
						//Print hash
						// printf("%lu\n", buf[j].x >> 8);
						//Print position
						// printf("%lu\n", ((buf[j].y << 32) >> 33) + 1 - k);

						int zp = (int) ((buf[j].y << 63) >> 63);
						buf[j].y = (buf[j].y >> 32) << 32 | (uint32_t)hcnt<<1 | zp;
						kv_push(mm128_t, km, *p, buf[j]);
						++hcnt;
					}
				}
				for (j = 0; j <= buf_pos; ++j){
					if (min.x == buf[j].x && min.y != buf[j].y && !bLstmers.contains(buf[j].x >> 8)){
						//Testing
						// if((buf[j].x >> 8) == 1389627116){
						// 	printf("i:%d\n", i);
						// 	// printf("Sixth push\n");
						// }
						// if(i >= 62396600 && i <=62396619){
						// printf("Hash in iteration i=%d: %lu\n", i, buf[j].x >> 8);
						// }
						//Print hash
						// printf("%lu\n", buf[j].x >> 8);
						//Print position
						// printf("%lu\n", ((buf[j].y << 32) >> 33) + 1 - k);

						int zp = (int) ((buf[j].y << 63) >> 63);
						buf[j].y = (buf[j].y >> 32) << 32 | (uint32_t)hcnt<<1 | zp;
						kv_push(mm128_t, km, *p, buf[j]);
						++hcnt;
					}
				}
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX && !bLstmers.contains(min.x >> 8)){
		//Testing
		// if((min.x >> 8) == 1389627116){
		// 	printf("i:%d\n", i);
		// 	// printf("Seventh push\n");
		// }
		// if(i >= 62396600 && i <=62396619){
		// printf("Hash in iteration i=%d: %lu\n", i, min.x >> 8);
		// }
		//Print hash
		// printf("%lu\n", min.x >> 8);
		//Print position
		// printf("%lu\n", ((min.y << 32) >> 33) + 1 - k);

		int zp = (int) ((min.y << 63) >> 63);
		min.y = (min.y >> 32) << 32 | (uint32_t)hcnt<<1 | zp;
		kv_push(mm128_t, km, *p, min);
		++hcnt;
	}
}
