#include <gtest/gtest.h>

#include "IOtest.h"
#include "IO.cpp"

//Tests for function const bool prsArgs(int&, char**, string&, string&, Measure&)//
//	1. Minimum number of parameters is (not) given DONE
//	2. First input sequence is (not) given DONE
//	3. Second input sequence is (not) given DONE
//	4. Similarity measure is (not) given DONE
//	5. Help flag is (not) set DONE
//	6. "intersec"/"algnWoutOffs" is given as a similarity measure DONE
//	7. "intersec" is (not) given first as a similarity measure and (not) again as a second DONE
//	8. "intersec" is (not) given first as a similarity measure, but (not) as a second DONE

//Tests the function prsArgs under the following conditions
//	1. Minimum number of parameters is not given
TEST_F(PrsArgsTest, fewPrms){
	nbArgs = 5;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("CalcSim");
	argv[1] = strdup("-a");
	argv[2] = strdup("-b");
	argv[3] = strdup("ACGT");
	argv[4] = strdup("-i");

	EXPECT_FALSE(prsArgs(nbArgs, argv, sa, sb, m));
	EXPECT_EQ(nbArgs, 5);
	EXPECT_EQ(sa, "");
	EXPECT_EQ(sb, "");
	EXPECT_EQ(none, m);
}

//Tests the function prsArgs under the following conditions
//	1. Minimum number of parameters is given
//	2. First input sequence is given
//	3. Second input sequence is given
//	4. Similarity measure is given
//	5. Help flag is set
//	6. "intersec" is given as a similarity measure
TEST_F(PrsArgsTest, hlpFlg){
	nbArgs = 7;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("CalcSim");
	argv[1] = strdup("-a");
	argv[2] = strdup("ACGT");
	argv[3] = strdup("-b");
	argv[4] = strdup("ACGT");
	argv[5] = strdup("-i");
	argv[6] = strdup("-h");

	EXPECT_FALSE(prsArgs(nbArgs, argv, sa, sb, m));
	EXPECT_EQ(nbArgs, 7);
	EXPECT_EQ(sa, "ACGT");
	EXPECT_EQ(sb, "ACGT");
	EXPECT_EQ(intersec, m);
}

//Tests the function prsArgs under the following conditions
//	1. Minimum number of parameters is given
//	2. First input sequence is not given
//	3. Second input sequence is given
//	4. Similarity measure is given
//	5. Help flag is not set
//	6. "intersec" is given as a similarity measure
TEST_F(PrsArgsTest, noSeqA){
	nbArgs = 10;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("CalcSim");
	argv[7] = strdup("-b");
	argv[8] = strdup("ACGT");
	argv[9] = strdup("-i");

	EXPECT_FALSE(prsArgs(nbArgs, argv, sa, sb, m));
	EXPECT_EQ(nbArgs, 10);
	EXPECT_EQ(sa, "");
	EXPECT_EQ(sb, "ACGT");
	EXPECT_EQ(intersec, m);
}

//Tests the function prsArgs under the following conditions
//	1. Minimum number of parameters is given
//	2. First input sequence is given
//	3. Second input sequence is not given
//	4. Similarity measure is given
//	5. Help flag is not set
//	6. "intersec" is given as a similarity measure
TEST_F(PrsArgsTest, noSeqB){
	nbArgs = 13;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("CalcSim");
	argv[10] = strdup("-a");
	argv[11] = strdup("ACGT");
	argv[12] = strdup("-i");

	EXPECT_FALSE(prsArgs(nbArgs, argv, sa, sb, m));
	EXPECT_EQ(nbArgs, 13);
	EXPECT_EQ(sa, "ACGT");
	EXPECT_EQ(sb, "");
	EXPECT_EQ(intersec, m);
}

//Tests the function prsArgs under the following conditions
//	1. Minimum number of parameters is given
//	2. First input sequence is given
//	3. Second input sequence is given
//	4. Similarity measure is not given
//	5. Help flag is not set
TEST_F(PrsArgsTest, noMes){
	nbArgs = 17;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("CalcSim");
	argv[13] = strdup("-b");
	argv[14] = strdup("ACGT");
	argv[15] = strdup("-a");
	argv[16] = strdup("ACGT");

	EXPECT_FALSE(prsArgs(nbArgs, argv, sa, sb, m));
	EXPECT_EQ(nbArgs, 17);
	EXPECT_EQ(sa, "ACGT");
	EXPECT_EQ(sb, "ACGT");
	EXPECT_EQ(none, m);
}

//Tests the function prsArgs under the following conditions
//	1. Minimum number of parameters is given
//	2. First input sequence is given
//	3. Second input sequence is given
//	4. Similarity measure is given
//	5. Help flag is not set
//	6. "algnWoutOffs" is given as a similarity measure
//	7. "intersec" is not given first as a similarity measure and not again as a second
TEST_F(PrsArgsTest, algMes){
	nbArgs = 23;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("CalcSim");
	argv[17] = strdup("-l");
	argv[18] = strdup("-a");
	argv[19] = strdup("ACGT");
	argv[20] = strdup("-b");
	argv[21] = strdup("ACGT");
	argv[22] = strdup("-l");

	EXPECT_TRUE(prsArgs(nbArgs, argv, sa, sb, m));
	EXPECT_EQ(nbArgs, 23);
	EXPECT_EQ(sa, "ACGT");
	EXPECT_EQ(sb, "ACGT");
	EXPECT_EQ(algnWoutOffs, m);
}

//Tests the function prsArgs under the following conditions
//	1. Minimum number of parameters is given
//	2. First input sequence is given
//	3. Second input sequence is given
//	4. Similarity measure is given
//	5. Help flag is not set
//	6. "intersec" is given as a similarity measure
//	7. "intersec" is given first as a similarity measure and again as a second
TEST_F(PrsArgsTest, dupInt){
	nbArgs = 29;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("CalcSim");
	argv[23] = strdup("-i");
	argv[24] = strdup("-a");
	argv[25] = strdup("ACGT");
	argv[26] = strdup("-b");
	argv[27] = strdup("ACGT");
	argv[28] = strdup("-i");

	EXPECT_TRUE(prsArgs(nbArgs, argv, sa, sb, m));
	EXPECT_EQ(nbArgs, 29);
	EXPECT_EQ(sa, "ACGT");
	EXPECT_EQ(sb, "ACGT");
	EXPECT_EQ(intersec, m);
}

//Tests the function prsArgs under the following conditions
//	1. Minimum number of parameters is given
//	2. First input sequence is given
//	3. Second input sequence is given
//	4. Similarity measure is given
//	5. Help flag is not set
//	6. "intersec"/"algnWoutOffs" is given as a similarity measure
//	8. "intersec" is given first as a similarity measure, but not as a second
TEST_F(PrsArgsTest, dupMes){
	nbArgs = 35;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("CalcSim");
	argv[29] = strdup("-i");
	argv[30] = strdup("-a");
	argv[31] = strdup("ACGT");
	argv[32] = strdup("-b");
	argv[33] = strdup("ACGT");
	argv[34] = strdup("-l");

	EXPECT_FALSE(prsArgs(nbArgs, argv, sa, sb, m));
	EXPECT_EQ(nbArgs, 35);
	EXPECT_EQ(sa, "ACGT");
	EXPECT_EQ(sb, "ACGT");
	EXPECT_EQ(intersec, m);
}

//Tests the function prsArgs under the following conditions
//	1. Minimum number of parameters is given
//	2. First input sequence is given
//	3. Second input sequence is given
//	4. Similarity measure is given
//	5. Help flag is not set
//	6. "intersec"/"algnWoutOffs" is given as a similarity measure
//	8. "intersec" is not given first as a similarity measure, but as a second
TEST_F(PrsArgsTest, dupOmes){
	nbArgs = 41;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[0] = strdup("CalcSim");
	argv[35] = strdup("-l");
	argv[36] = strdup("-a");
	argv[37] = strdup("ACGT");
	argv[38] = strdup("-b");
	argv[39] = strdup("ACGT");
	argv[40] = strdup("-i");

	EXPECT_FALSE(prsArgs(nbArgs, argv, sa, sb, m));
	EXPECT_EQ(nbArgs, 41);
	EXPECT_EQ(sa, "ACGT");
	EXPECT_EQ(sb, "ACGT");
	EXPECT_EQ(algnWoutOffs, m);
}
