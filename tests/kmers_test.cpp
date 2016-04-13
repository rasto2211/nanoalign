#include <vector>
#include <string>

#include "src/kmers.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using testing::StrEq;

TEST(KmersTest, kmersInDist0Test) {
  std::vector<std::string> res = kmersInDist("ACTGC", 0);
  EXPECT_THAT(res, ::testing::ElementsAreArray({StrEq("ACTGC")}));
}

TEST(KmersTest, kmersInDist1Test) {
  std::vector<std::string> res = kmersInDist("ACTGC", 1);
  EXPECT_THAT(res,
              ::testing::ElementsAreArray({StrEq("CTGCA"), StrEq("CTGCC"),
                                           StrEq("CTGCT"), StrEq("CTGCG")}));
}

TEST(KmersTest, kmersInDist2Test) {
  std::vector<std::string> res = kmersInDist("ACTGC", 2);
  EXPECT_THAT(
      res,
      ::testing::ElementsAreArray(
          {StrEq("TGCAA"), StrEq("TGCAC"), StrEq("TGCAT"), StrEq("TGCAG"),
           StrEq("TGCCA"), StrEq("TGCCC"), StrEq("TGCCT"), StrEq("TGCCG"),
           StrEq("TGCTA"), StrEq("TGCTC"), StrEq("TGCTT"), StrEq("TGCTG"),
           StrEq("TGCGA"), StrEq("TGCGC"), StrEq("TGCGT"), StrEq("TGCGG")}));
}

TEST(KmersTest, kmersInDistTooFarTes) {
  std::vector<std::string> res = kmersInDist("AC", 4);
  EXPECT_THAT(res, ::testing::ElementsAreArray(
                       {StrEq("AA"), StrEq("AC"), StrEq("AT"), StrEq("AG"),
                        StrEq("CA"), StrEq("CC"), StrEq("CT"), StrEq("CG"),
                        StrEq("TA"), StrEq("TC"), StrEq("TT"), StrEq("TG"),
                        StrEq("GA"), StrEq("GC"), StrEq("GT"), StrEq("GG")}));
}

TEST(KmersTest, KmerStrToIntLongLongTest) {
  long long num_kmer = encodeKmer<long long>("TTCGGTTCGACGTTGACCTCCATTATCT");
  EXPECT_EQ(119320958947904038LL, num_kmer);
}

TEST(KmersTest, KmerStrToIntTest) {
  int num_kmer = encodeKmer<int>("ACTG");
  EXPECT_EQ(283, num_kmer);
}

TEST(KmersTest, KmerStrToInt14LengthTest) {
  int num_kmer = encodeKmer<int>("GGGGGGGGGGGGGG");
  EXPECT_EQ(536870911, num_kmer);
}

TEST(KmersTest, IntToKmerTest) {
  std::string kmer = decodeKmer<int>(283);
  EXPECT_EQ("ACTG", kmer);
}

TEST(KmersTest, LongLongToKmerTest) {
  std::string kmer = decodeKmer<long long>(119320958947904038LL);
  EXPECT_EQ("TTCGGTTCGACGTTGACCTCCATTATCT", kmer);
}

TEST(KmersTest, WindowIteratorTest) {
  std::string input_seq = "AACTGATC";
  KmerWindowIterator<int> window_it =
      KmerWindowIterator<int>(5, input_seq.begin(), input_seq.end());

  // Initial window.
  EXPECT_EQ(encodeKmer<int>("AACTG"), window_it.currentKmerCode());
  EXPECT_EQ("AACTG", window_it.currentKmer());
  EXPECT_TRUE(window_it.hasNext());

  // Next window.
  EXPECT_EQ(encodeKmer<int>("ACTGA"), window_it.next());
  EXPECT_EQ("ACTGA", window_it.currentKmer());
  EXPECT_TRUE(window_it.hasNext());

  // Next window.
  EXPECT_EQ(encodeKmer<int>("CTGAT"), window_it.next());
  EXPECT_EQ("CTGAT", window_it.currentKmer());
  EXPECT_TRUE(window_it.hasNext());

  // Next window.
  EXPECT_EQ(encodeKmer<int>("TGATC"), window_it.next());
  EXPECT_EQ("TGATC", window_it.currentKmer());
  EXPECT_FALSE(window_it.hasNext());

  EXPECT_EQ(-1, window_it.next());
}

TEST(KmersTest, NumKmersOfLengthTest) {
  EXPECT_EQ(1024, numKmersOf(5));  // 4^5
}

TEST(KmersTest, KmerToLexicographicPosTest) {
  EXPECT_EQ(1, kmerToLexicographicPos("AAA"));
  EXPECT_EQ(2, kmerToLexicographicPos("AAC"));
  EXPECT_EQ(3, kmerToLexicographicPos("AAT"));
  EXPECT_EQ(4, kmerToLexicographicPos("AAG"));
  EXPECT_EQ(5, kmerToLexicographicPos("ACA"));
  EXPECT_EQ(6, kmerToLexicographicPos("ACC"));
  EXPECT_EQ(7, kmerToLexicographicPos("ACT"));
  EXPECT_EQ(8, kmerToLexicographicPos("ACG"));

  EXPECT_EQ(64, kmerToLexicographicPos("GGG"));
}

TEST(KmersTest, KmerInLexicographicPosTest) {
  EXPECT_EQ("AAA", kmerInLexicographicPos(1, 3));
  EXPECT_EQ("AAC", kmerInLexicographicPos(2, 3));
  EXPECT_EQ("AAT", kmerInLexicographicPos(3, 3));
  EXPECT_EQ("AAG", kmerInLexicographicPos(4, 3));
  EXPECT_EQ("ACA", kmerInLexicographicPos(5, 3));
  EXPECT_EQ("ACC", kmerInLexicographicPos(6, 3));
  EXPECT_EQ("ACT", kmerInLexicographicPos(7, 3));
  EXPECT_EQ("ACG", kmerInLexicographicPos(8, 3));

  EXPECT_EQ("GGG", kmerInLexicographicPos(64, 3));
}

TEST(KmersTest, kmersUpToDist2Test) {
  std::unordered_set<std::string> res = kmersUpToDist("ACTGC", 2);
  EXPECT_THAT(
      res, ::testing::UnorderedElementsAreArray(
               {StrEq("ACTGC"), StrEq("CTGCA"), StrEq("CTGCC"), StrEq("CTGCT"),
                StrEq("CTGCG"), StrEq("TGCAA"), StrEq("TGCAC"), StrEq("TGCAT"),
                StrEq("TGCAG"), StrEq("TGCCA"), StrEq("TGCCC"), StrEq("TGCCT"),
                StrEq("TGCCG"), StrEq("TGCTA"), StrEq("TGCTC"), StrEq("TGCTT"),
                StrEq("TGCTG"), StrEq("TGCGA"), StrEq("TGCGC"), StrEq("TGCGT"),
                StrEq("TGCGG")}));
}

TEST(KmersTest, kmersUpToDistHomopolymerTest) {
  std::unordered_set<std::string> res = kmersUpToDist("AAA", 2);
  EXPECT_THAT(res,
              ::testing::UnorderedElementsAreArray(
                  {"AAA", "AAC", "AAT", "AAG", "ACA", "ACT", "ACC", "ACG",
                   "ATA", "ATC", "ATT", "ATG", "AGA", "AGC", "AGT", "AGG"}));
}

TEST(KmersTest, GetMoveTest) {
  EXPECT_EQ(-1, getMove("GTTCG", "GTTC"));
  EXPECT_EQ(0, getMove("GTTCG", "GTTCG"));
  EXPECT_EQ(1, getMove("GTTCG", "TTCGG"));
  EXPECT_EQ(2, getMove("GTTCG", "TCGGA"));
  EXPECT_EQ(3, getMove("GTTCG", "CGGAA"));
  EXPECT_EQ(4, getMove("GTTCG", "GGAAT"));
  EXPECT_EQ(5, getMove("GTTCG", "AATTC"));
}
