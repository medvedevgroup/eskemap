#include <gtest/gtest.h>

#include "IOtest.h"

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

//Tests for function const bool prsArgs(int&, char**, string&, string&, uint32_t&, uint32_t&, double&, string&, uint32_t&, float&, 
//float&, bool&, float&, float&, bool&)//
//	1. Pattern sequence is (not) given DONE
//	2. Text sequence is (not) given DONE
//	3. K-mer length is (not) given DONE
//	4. K-mer length is given and is (not) positive DONE
//	5. Hash ratio is (not) given DONE
//	6. Hash ratio is given and is (not) positive DONE
//	7. Hash ratio is given and is (not) larger than MAX_RATIO DONE
//	8. Common hash weight is (not) given DONE
//	9. Common hash weight is given and is (not) positive DONE
//	10. Unique hash weight is (not) given DONE
//	11. Unique hash weight is given and is (not) positive DONE
//	12. t-homology threshold is (not) given DONE
//	13. Normalization flag is (not) given DONE
//	14. Help flag is (not) given DONE
//	15. Blacklist is (not) given DONE

//Tests the function prsArgs under the following conditions
//	1. Pattern sequence is given
//	2. Text sequence is given
//	3. K-mer length is given
//	4. K-mer length is given and is positive
//	5. Hash ratio is given
//	6. Hash ratio is given and is positive
//	7. Hash ratio is given and is not larger than MAX_RATIO
//	8. Common hash weight is given
//	9. Common hash weight is given and is positive
//	10. Unique hash weight is given
//	11. Unique hash weight is given and is positive
//	12. t-homology threshold is given
//	13. Normalization flag is given
//	14. Help flag is not given
//	15. Blacklist is not given
TEST_F(PrsArgs1Test, pttnGvn){
	nbArgs = 57;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[41] = strdup("FindThoms");
	argv[42] = strdup("-p");
	argv[43] = strdup("A");
	argv[44] = strdup("-s");
	argv[45] = strdup("ACGT");
	argv[46] = strdup("-r");
	argv[47] = strdup("0.2");
	argv[48] = strdup("-c");
	argv[49] = strdup("2");
	argv[50] = strdup("-u");
	argv[51] = strdup("2");
	argv[52] = strdup("-t");
	argv[53] = strdup("-1");
	argv[54] = strdup("-n");
	argv[55] = strdup("-k");
	argv[56] = strdup("2147483647");

	EXPECT_TRUE(prsArgs(nbArgs, argv, p, s, k, h, b, c, u, t, n));
	EXPECT_EQ(nbArgs, 57);
	EXPECT_EQ(p, "A");
	EXPECT_EQ(s, "ACGT");
	EXPECT_EQ(k, 2147483647);
	EXPECT_EQ(h, 0.2);
	EXPECT_EQ(c, 2);
	EXPECT_EQ(u, 2);
	EXPECT_EQ(t, -1);
	EXPECT_TRUE(n);
	EXPECT_TRUE(b.empty());
}

//Tests the function prsArgs under the following conditions
//	1. Pattern sequence is not given
//	2. Text sequence is given
//	3. K-mer length is not given
//	5. Hash ratio is not given
//	8. Common hash weight is not given
//	10. Unique hash weight is not given
//	12. t-homology threshold is not given
//	13. Normalization flag is not given
//	14. Help flag is not given
//	15. Blacklist is not given
TEST_F(PrsArgs1Test, noPttn){
	nbArgs = 59;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[56] = strdup("FindThoms");
	argv[57] = strdup("-s");
	argv[58] = strdup("A");

	EXPECT_FALSE(prsArgs(nbArgs, argv, p, s, k, h, b, c, u, t, n));
	EXPECT_EQ(nbArgs, 59);
	EXPECT_EQ(p, "");
	EXPECT_EQ(s, "A");
	EXPECT_EQ(k, K);
	EXPECT_EQ(h, HASH_RATIO);
	EXPECT_EQ(c, DEFAULT_WEIGHT);
	EXPECT_EQ(u, DEFAULT_WEIGHT);
	EXPECT_EQ(t, T);
	EXPECT_FALSE(n);
	EXPECT_TRUE(b.empty());
}

//Tests the function prsArgs under the following conditions
//	1. Pattern sequence is given
//	2. Text sequence is not given
//	3. K-mer length is not given
//	5. Hash ratio is not given
//	8. Common hash weight is not given
//	10. Unique hash weight is not given
//	12. t-homology threshold is not given
//	13. Normalization flag is not given
//	14. Help flag is not given
//	15. Blacklist is not given
TEST_F(PrsArgs1Test, noTxt){
	nbArgs = 61;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[58] = strdup("FindThoms");
	argv[59] = strdup("-p");
	argv[60] = strdup("A");

	EXPECT_FALSE(prsArgs(nbArgs, argv, p, s, k, h, b, c, u, t, n));
	EXPECT_EQ(nbArgs, 61);
	EXPECT_EQ(p, "A");
	EXPECT_EQ(s, "");
	EXPECT_EQ(k, K);
	EXPECT_EQ(h, HASH_RATIO);
	EXPECT_EQ(c, DEFAULT_WEIGHT);
	EXPECT_EQ(u, DEFAULT_WEIGHT);
	EXPECT_EQ(t, T);
	EXPECT_FALSE(n);
	EXPECT_TRUE(b.empty());
}

//Tests the function prsArgs under the following conditions
//	1. Pattern sequence is not given
//	2. Text sequence is not given
//	3. K-mer length is not given
//	5. Hash ratio is not given
//	8. Common hash weight is not given
//	10. Unique hash weight is not given
//	12. t-homology threshold is not given
//	13. Normalization flag is not given
//	14. Help flag is not given
//	15. Blacklist is not given
TEST_F(PrsArgs1Test, noSeqs){
	nbArgs = 61;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[60] = strdup("FindThoms");

	EXPECT_FALSE(prsArgs(nbArgs, argv, p, s, k, h, b, c, u, t, n));
	EXPECT_EQ(nbArgs, 61);
	EXPECT_EQ(p, "");
	EXPECT_EQ(s, "");
	EXPECT_EQ(k, K);
	EXPECT_EQ(h, HASH_RATIO);
	EXPECT_EQ(c, DEFAULT_WEIGHT);
	EXPECT_EQ(u, DEFAULT_WEIGHT);
	EXPECT_EQ(t, T);
	EXPECT_FALSE(n);
	EXPECT_TRUE(b.empty());
}

//Tests the function prsArgs under the following conditions
//	1. Pattern sequence is not given
//	2. Text sequence is not given
//	3. K-mer length is given
//	4. K-mer length is given and is not positive
//	5. Hash ratio is not given
//	8. Common hash weight is not given
//	10. Unique hash weight is not given
//	12. t-homology threshold is not given
//	13. Normalization flag is not given
//	14. Help flag is not given
//	15. Blacklist is not given
TEST_F(PrsArgs1Test, negK){
	nbArgs = 63;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[60] = strdup("FindThoms");
	argv[61] = strdup("-k");
	argv[62] = strdup("-1");

	EXPECT_FALSE(prsArgs(nbArgs, argv, p, s, k, h, b, c, u, t, n));
	EXPECT_EQ(nbArgs, 63);
	EXPECT_EQ(p, "");
	EXPECT_EQ(s, "");
	EXPECT_EQ(k, K);
	EXPECT_EQ(h, HASH_RATIO);
	EXPECT_EQ(c, DEFAULT_WEIGHT);
	EXPECT_EQ(u, DEFAULT_WEIGHT);
	EXPECT_EQ(t, T);
	EXPECT_FALSE(n);
	EXPECT_TRUE(b.empty());
}

//Tests the function prsArgs under the following conditions
//	1. Pattern sequence is not given
//	2. Text sequence is not given
//	3. K-mer length is not given
//	5. Hash ratio is given
//	6. Hash ratio is given and is not positive
//	7. Hash ratio is given and is not larger than MAX_RATIO
//	8. Common hash weight is not given
//	10. Unique hash weight is not given
//	12. t-homology threshold is not given
//	13. Normalization flag is not given
//	14. Help flag is not given
//	15. Blacklist is not given
TEST_F(PrsArgs1Test, negR){
	nbArgs = 65;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[62] = strdup("FindThoms");
	argv[63] = strdup("-r");
	argv[64] = strdup("-1");

	EXPECT_FALSE(prsArgs(nbArgs, argv, p, s, k, h, b, c, u, t, n));
	EXPECT_EQ(nbArgs, 65);
	EXPECT_EQ(p, "");
	EXPECT_EQ(s, "");
	EXPECT_EQ(k, K);
	EXPECT_EQ(h, HASH_RATIO);
	EXPECT_EQ(c, DEFAULT_WEIGHT);
	EXPECT_EQ(u, DEFAULT_WEIGHT);
	EXPECT_EQ(t, T);
	EXPECT_FALSE(n);
	EXPECT_TRUE(b.empty());
}

//Tests the function prsArgs under the following conditions
//	1. Pattern sequence is not given
//	2. Text sequence is not given
//	3. K-mer length is not given
//	5. Hash ratio is given
//	6. Hash ratio is given and is positive
//	7. Hash ratio is given and is larger than MAX_RATIO
//	8. Common hash weight is not given
//	10. Unique hash weight is not given
//	12. t-homology threshold is not given
//	13. Normalization flag is not given
//	14. Help flag is not given
//	15. Blacklist is not given
TEST_F(PrsArgs1Test, lrgR){
	nbArgs = 67;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[64] = strdup("FindThoms");
	argv[65] = strdup("-r");
	argv[66] = strdup("1.1");

	EXPECT_FALSE(prsArgs(nbArgs, argv, p, s, k, h, b, c, u, t, n));
	EXPECT_EQ(nbArgs, 67);
	EXPECT_EQ(p, "");
	EXPECT_EQ(s, "");
	EXPECT_EQ(k, K);
	EXPECT_EQ(h, HASH_RATIO);
	EXPECT_EQ(c, DEFAULT_WEIGHT);
	EXPECT_EQ(u, DEFAULT_WEIGHT);
	EXPECT_EQ(t, T);
	EXPECT_FALSE(n);
	EXPECT_TRUE(b.empty());
}

//Tests the function prsArgs under the following conditions
//	1. Pattern sequence is not given
//	2. Text sequence is not given
//	3. K-mer length is not given
//	5. Hash ratio is not given
//	8. Common hash weight is given
//	9. Common hash weight is given and is not positive
//	10. Unique hash weight is not given
//	12. t-homology threshold is not given
//	13. Normalization flag is not given
//	14. Help flag is not given
//	15. Blacklist is not given
TEST_F(PrsArgs1Test, negC){
	nbArgs = 69;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[66] = strdup("FindThoms");
	argv[67] = strdup("-c");
	argv[68] = strdup("-1");

	EXPECT_FALSE(prsArgs(nbArgs, argv, p, s, k, h, b, c, u, t, n));
	EXPECT_EQ(nbArgs, 69);
	EXPECT_EQ(p, "");
	EXPECT_EQ(s, "");
	EXPECT_EQ(k, K);
	EXPECT_EQ(h, HASH_RATIO);
	EXPECT_EQ(c, DEFAULT_WEIGHT);
	EXPECT_EQ(u, DEFAULT_WEIGHT);
	EXPECT_EQ(t, T);
	EXPECT_FALSE(n);
	EXPECT_TRUE(b.empty());
}

//Tests the function prsArgs under the following conditions
//	1. Pattern sequence is not given
//	2. Text sequence is not given
//	3. K-mer length is not given
//	5. Hash ratio is not given
//	8. Common hash weight is not given
//	10. Unique hash weight is given
//	11. Unique hash weight is given and is not positive
//	12. t-homology threshold is not given
//	13. Normalization flag is not given
//	14. Help flag is not given
//	15. Blacklist is not given
TEST_F(PrsArgs1Test, negU){
	nbArgs = 71;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[68] = strdup("FindThoms");
	argv[69] = strdup("-u");
	argv[70] = strdup("-1");

	EXPECT_FALSE(prsArgs(nbArgs, argv, p, s, k, h, b, c, u, t, n));
	EXPECT_EQ(nbArgs, 71);
	EXPECT_EQ(p, "");
	EXPECT_EQ(s, "");
	EXPECT_EQ(k, K);
	EXPECT_EQ(h, HASH_RATIO);
	EXPECT_EQ(c, DEFAULT_WEIGHT);
	EXPECT_EQ(u, DEFAULT_WEIGHT);
	EXPECT_EQ(t, T);
	EXPECT_FALSE(n);
	EXPECT_TRUE(b.empty());
}

//Tests the function prsArgs under the following conditions
//	1. Pattern sequence is not given
//	2. Text sequence is not given
//	3. K-mer length is not given
//	5. Hash ratio is not given
//	8. Common hash weight is not given
//	10. Unique hash weight is not given
//	12. t-homology threshold is not given
//	13. Normalization flag is not given
//	14. Help flag is given
//	15. Blacklist is not given
TEST_F(PrsArgs1Test, help){
	nbArgs = 72;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[70] = strdup("FindThoms");
	argv[71] = strdup("-h");

	EXPECT_FALSE(prsArgs(nbArgs, argv, p, s, k, h, b, c, u, t, n));
	EXPECT_EQ(nbArgs, 72);
	EXPECT_EQ(p, "");
	EXPECT_EQ(s, "");
	EXPECT_EQ(k, K);
	EXPECT_EQ(h, HASH_RATIO);
	EXPECT_EQ(c, DEFAULT_WEIGHT);
	EXPECT_EQ(u, DEFAULT_WEIGHT);
	EXPECT_EQ(t, T);
	EXPECT_FALSE(n);
	EXPECT_TRUE(b.empty());
}

//Tests the function prsArgs under the following conditions
//	1. Pattern sequence is not given
//	2. Text sequence is not given
//	3. K-mer length is not given
//	5. Hash ratio is not given
//	8. Common hash weight is not given
//	10. Unique hash weight is not given
//	12. t-homology threshold is not given
//	13. Normalization flag is not given
//	14. Help flag is not given
//	15. Blacklist is given
TEST_F(PrsArgs1Test, blGvn){
	nbArgs = 74;
	argv = (char**) malloc(nbArgs * sizeof(char*));
	argv[71] = strdup("FindThoms");
	argv[72] = strdup("-b");
	argv[73] = strdup("testBlacklist.txt");

	EXPECT_FALSE(prsArgs(nbArgs, argv, p, s, k, h, b, c, u, t, n));
	EXPECT_EQ(nbArgs, 74);
	EXPECT_EQ(p, "");
	EXPECT_EQ(s, "");
	EXPECT_EQ(k, K);
	EXPECT_EQ(h, HASH_RATIO);
	EXPECT_EQ(c, DEFAULT_WEIGHT);
	EXPECT_EQ(u, DEFAULT_WEIGHT);
	EXPECT_EQ(t, T);
	EXPECT_FALSE(n);
	EXPECT_EQ(b, "testBlacklist.txt");
}

//Tests for function const bool readFASTA(const string&, string&)//
//	1. We can(not) open the file DONE
//	2. We have (not) found a second header line DONE
//	3. We have (not) found the beginning of a header line DONE
//	4. We have (not) found a line break DONE
//	5. We are (not) inside a header line DONE
//	6. We have (not) read a valid DNA nucleotide DONE
//	7. We have found a '>' inside a header DONE

//Tests the function readFASTA under the following conditions
//	1. We can open the file
//	2. We have (not) found a second header line
//	3. We have (not) found the beginning of a header line
//	4. We have (not) found a line break
//	5. We are (not) inside a header line
//	6. We have (not) read a valid DNA nucleotide
//	7. We have found a '>' inside a header
TEST_F(ReadFASTAtest, valFl){
	EXPECT_TRUE(readFASTA("TestFASTA.fasta", s));
	EXPECT_EQ(s, "AGT");
}

//Tests the function readFASTA under the following conditions
//	1. We cannot open the file
TEST_F(ReadFASTAtest, invFl){
	EXPECT_FALSE(readFASTA("A", s));
	EXPECT_EQ(s, "");
}

//Tests for function inline void outputHoms(const vector<Thomology>&, const bool&, const uint32_t&)//
//	1. Scores are (not) normalized DONE
//	2. Scores are normalized and substring is (not) longer than pattern DONE

//Tests the function outputHoms under the following conditions
//	1. Scores are normalized
//	2. Scores are normalized and substring is (not) longer than pattern
TEST_F(OutputHomsTest, nScrs){
	r.push_back(make_tuple(1, 4, 2));
	r.push_back(make_tuple(1, 2, 2));
	r.push_back(make_tuple(7, 9, 3));

	coutPtr = cout.rdbuf();
	cout.rdbuf(stream.rdbuf());

	outputHoms(r, true, 3);

	stream >> res;
	EXPECT_EQ("i:", res);
	stream >> res;
	EXPECT_EQ("1", res);
	stream >> res;
	EXPECT_EQ("j:", res);
	stream >> res;
	EXPECT_EQ("4", res);
	stream >> res;
	EXPECT_EQ("score:", res);
	stream >> res;
	EXPECT_EQ("0.5", res);
	stream >> res;
	EXPECT_EQ("i:", res);
	stream >> res;
	EXPECT_EQ("1", res);
	stream >> res;
	EXPECT_EQ("j:", res);
	stream >> res;
	EXPECT_EQ("2", res);
	stream >> res;
	EXPECT_EQ("score:", res);
	stream >> res;
	EXPECT_EQ("0.666667", res);
	stream >> res;
	EXPECT_EQ("i:", res);
	stream >> res;
	EXPECT_EQ("7", res);
	stream >> res;
	EXPECT_EQ("j:", res);
	stream >> res;
	EXPECT_EQ("9", res);
	stream >> res;
	EXPECT_EQ("score:", res);
	stream >> res;
	EXPECT_EQ("1", res);
	stream >> res;
	EXPECT_EQ("1", res);
	EXPECT_TRUE(stream.eof());

	cout.rdbuf(coutPtr);
}

//Tests the function outputHoms under the following conditions
//	1. Scores are not normalized
TEST_F(OutputHomsTest, nnScrs){
	r.push_back(make_tuple(1, 4, 2));
	r.push_back(make_tuple(1, 2, 2));
	r.push_back(make_tuple(7, 9, 3));

	coutPtr = cout.rdbuf();
	cout.rdbuf(stream.rdbuf());

	outputHoms(r, false, 3);

	stream >> res;
	EXPECT_EQ("i:", res);
	stream >> res;
	EXPECT_EQ("1", res);
	stream >> res;
	EXPECT_EQ("j:", res);
	stream >> res;
	EXPECT_EQ("4", res);
	stream >> res;
	EXPECT_EQ("score:", res);
	stream >> res;
	EXPECT_EQ("2", res);
	stream >> res;
	EXPECT_EQ("i:", res);
	stream >> res;
	EXPECT_EQ("1", res);
	stream >> res;
	EXPECT_EQ("j:", res);
	stream >> res;
	EXPECT_EQ("2", res);
	stream >> res;
	EXPECT_EQ("score:", res);
	stream >> res;
	EXPECT_EQ("2", res);
	stream >> res;
	EXPECT_EQ("i:", res);
	stream >> res;
	EXPECT_EQ("7", res);
	stream >> res;
	EXPECT_EQ("j:", res);
	stream >> res;
	EXPECT_EQ("9", res);
	stream >> res;
	EXPECT_EQ("score:", res);
	stream >> res;
	EXPECT_EQ("3", res);
	stream >> res;
	EXPECT_EQ("3", res);
	EXPECT_TRUE(stream.eof());

	cout.rdbuf(coutPtr);
}

//Tests for function const bool lPttnSks(ifstream&, const uint32_t&, const double&, vector<pair<string, Sketch>>&)//
//	1. The file we want to read is (not) open DONE
//	2. We have (not) found the beginning of a second sequence entry DONE
//	3. We found another sequence entry and did (not) reach the batch size limit DONE
//	4. We did (not) find a character to update the current sequence id DONE
//	5. We do (not) set the "header read" flag, because it was set before DONE
//	6. We do (not) set the "header read" flag, because we discover the beginning of a sequence entry for the first time DONE
//	7. We do (not) set the "id read" flag, because it was set before DONE
//	8. We do (not) set the "id read" flag and it was not set before DONE
//	9. We do (not) set the "line break discovered" flag, because it was set before DONE
//	10. We do (not) set the "line break discovered" flag and is was not set before DONE
//	11. We are (not) inside the header of a sequence entry DONE
//	12. We do (not) read another nucleotide DONE
//	13. We could (not) successfully generate a sketch of the last sequence entry's sequence DONE

//Tests the function lPttnSks under the following conditions
//	1. The file we want to read is open
//	2. We have (not) found the beginning of a second sequence entry
//	3. We found another sequence entry and did not reach the batch size limit
//	4. We did (not) find a character to update the current sequence id
//	5. We do (not) set the "header read" flag, because it was set before
//	6. We do (not) set the "header read" flag, because we discover the beginning of a sequence entry for the first time
//	7. We do (not) set the "id read" flag, because it was set before
//	8. We do (not) set the "id read" flag and it was not set before
//	9. We do (not) set the "line break discovered" flag, because it was set before
//	10. We do (not) set the "line break discovered" flag and is was not set before
//	11. We are (not) inside the header of a sequence entry
//	12. We do (not) read another nucleotide
//	13. We could not successfully generate a sketch of the last sequence entry's sequence
TEST_F(lPttnSksTest, exFl){
	fStr.open("TestFASTA.fasta");

	EXPECT_FALSE(lPttnSks(fStr, K, HASH_RATIO, b, s));

	ASSERT_FALSE(s.empty());
	EXPECT_EQ(s.size(), 1);
	EXPECT_EQ(s.front().first, "Header");
	EXPECT_TRUE(s.front().second.empty());
}

//Tests the function lPttnSks under the following conditions
//	1. The file we want to read is not open
TEST_F(lPttnSksTest, nexFl){
	fStr.open("TestFASTA.fsta");

	EXPECT_FALSE(lPttnSks(fStr, K, HASH_RATIO, b, s));

	EXPECT_TRUE(s.empty());
}

//Tests the function lPttnSks under the following conditions
//	1. The file we want to read is open
//	2. We have (not) found the beginning of a second sequence entry
//	3. We found another sequence entry and did (not) reach the batch size limit
//	4. We did (not) find a character to update the current sequence id
//	5. We do (not) set the "header read" flag, because it was set before
//	6. We do (not) set the "header read" flag, because we discover the beginning of a sequence entry for the first time
//	7. We do not set the "id read" flag, because it was set before
//	8. We do not set the "id read" flag and it was not set before
//	9. We do (not) set the "line break discovered" flag, because it was set before
//	10. We do (not) set the "line break discovered" flag and is was not set before
//	11. We are (not) inside the header of a sequence entry
//	12. We do (not) read another nucleotide
TEST_F(lPttnSksTest, btchLim){
	fStr.open("TestFASTA2.fasta");

	EXPECT_TRUE(lPttnSks(fStr, K, HASH_RATIO, b, s));

	ASSERT_FALSE(s.empty());
	EXPECT_EQ(s.size(), PATTERN_BATCH_SIZE);
	i = s.begin();

	for(uint32_t j = 0; j < PATTERN_BATCH_SIZE; ++j, ++i){
		EXPECT_EQ(i->first, "s" + to_string(j));
		EXPECT_TRUE(i->second.empty());
	}
}

//Tests the function lPttnSks under the following conditions
//	1. The file we want to read is open
//	2. We have (not) found the beginning of a second sequence entry
//	3. We found another sequence entry and did not reach the batch size limit
//	4. We did (not) find a character to update the current sequence id
//	5. We do (not) set the "header read" flag, because it was set before
//	6. We do (not) set the "header read" flag, because we discover the beginning of a sequence entry for the first time
//	7. We do not set the "id read" flag, because it was set before
//	8. We do not set the "id read" flag and it was not set before
//	9. We do (not) set the "line break discovered" flag, because it was set before
//	10. We do (not) set the "line break discovered" flag and is was not set before
//	11. We are (not) inside the header of a sequence entry
//	12. We do (not) read another nucleotide
//	13. We could successfully generate a sketch of the last sequence entry's sequence
TEST_F(lPttnSksTest, sucLseq){
	fStr.open("TestFASTA1.fasta");

	EXPECT_FALSE(lPttnSks(fStr, K, HASH_RATIO, b, s));

	ASSERT_FALSE(s.empty());
	i = s.begin();
	EXPECT_EQ(i->first, "testEntry1");
	EXPECT_EQ(i->second.size(), 1);
	EXPECT_EQ(i->second.front(), getHash(calcKmerNb("CGTACGTAC"), (2 << (2 * K) - 1) - 1));
	++i;
	EXPECT_EQ(i->first, "testEntry2");
	EXPECT_TRUE(i->second.empty());
	EXPECT_EQ(s.back().first, "testEntry3");
	EXPECT_TRUE(s.back().second.empty());
}

//Tests for function const unordered_map<uint64_t, char> readBlstKmers(const string&)//
//	1. The file we want to read is (not) open DONE
//	2. The file we want to read is (not) empty DONE

//Tests the function realBlstKmers under the following conditions
//	1. The file we want to read is open
//	2. The file we want to read is empty
TEST_F(ReadBlstKmersTest, empFl){
	b = readBlstKmers("emptyBlacklist.txt");

	EXPECT_TRUE(b.empty());
}

//Tests the function realBlstKmers under the following conditions
//	1. The file we want to read is not open
TEST_F(ReadBlstKmersTest, nExFl){
	b = readBlstKmers("emptyBlackList.txt");

	EXPECT_TRUE(b.empty());
}

//Tests the function realBlstKmers under the following conditions
//	1. The file we want to read is open
//	2. The file we want to read is not empty
TEST_F(ReadBlstKmersTest, nEmpLst){
	b = readBlstKmers("testBlacklist.txt");

	EXPECT_FALSE(b.empty());
	EXPECT_EQ(b.size(), 1);
	EXPECT_EQ(b[0], 1);
}
