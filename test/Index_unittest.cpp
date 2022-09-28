#include <gtest/gtest.h>

#include "IndexTest.h"
#include "Index.cpp"

//Tests for function inline const bool smHshnPos(const pair<uint64_t, uint32_t>&, const pair<uint64_t, uint32_t>&)//
//	1. Hashes are equal and first pair's position is (not) smaller 0/0
//	2. Hashes are not equal and first pair's position is (not) smaller 0/0

//Tests the function smHshnPos under the following conditions
//	1. Hashes are equal and first pair's position is smaller
TEST_F(IndexTest, hefs){
	fp = make_pair(42, 20);
	sp = make_pair(42, 21);

	EXPECT_TRUE(smHshnPos(fp, sp));
}
