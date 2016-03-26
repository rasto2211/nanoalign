#include <vector>
#include <string>

#include "kmers.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using testing::StrEq;

TEST(Kmers, AllNextKmersDist0Test) {
  std::vector<std::string> res = allNextKmers("ACTGC", 0);
  EXPECT_THAT(res, ::testing::ElementsAreArray({StrEq("ACTGC")}));
}

TEST(Kmers, AllNextKmersDist1Test) {
  std::vector<std::string> res = allNextKmers("ACTGC", 1);
  EXPECT_THAT(res,
              ::testing::ElementsAreArray({StrEq("CTGCA"), StrEq("CTGCC"),
                                           StrEq("CTGCT"), StrEq("CTGCG")}));
}

TEST(Kmers, AllNextKmersDist2Test) {
  std::vector<std::string> res = allNextKmers("ACTGC", 2);
  EXPECT_THAT(
      res,
      ::testing::ElementsAreArray(
          {StrEq("TGCAA"), StrEq("TGCAC"), StrEq("TGCAT"), StrEq("TGCAG"),
           StrEq("TGCCA"), StrEq("TGCCC"), StrEq("TGCCT"), StrEq("TGCCG"),
           StrEq("TGCTA"), StrEq("TGCTC"), StrEq("TGCTT"), StrEq("TGCTG"),
           StrEq("TGCGA"), StrEq("TGCGC"), StrEq("TGCGT"), StrEq("TGCGG")}));
}

TEST(Kmers, AllNextKmersDistTooMuchTest) {
  std::vector<std::string> res = allNextKmers("AC", 4);
  EXPECT_THAT(res, ::testing::ElementsAreArray(
                       {StrEq("AA"), StrEq("AC"), StrEq("AT"), StrEq("AG"),
                        StrEq("CA"), StrEq("CC"), StrEq("CT"), StrEq("CG"),
                        StrEq("TA"), StrEq("TC"), StrEq("TT"), StrEq("TG"),
                        StrEq("GA"), StrEq("GC"), StrEq("GT"), StrEq("GG")}));
}

TEST(Kmers, KmerStrToIntLongLongTest) {
  long long num_kmer = encodeKmer<long long>("TTCGGTTCGACGTTGACCTCCATTATCT");
  EXPECT_EQ(119320958947904038LL, num_kmer);
}

TEST(Kmers, KmerStrToIntTest) {
  int num_kmer = encodeKmer<int>("ACTG");
  EXPECT_EQ(283, num_kmer);
}

TEST(Kmers, KmerStrToInt14LengthTest) {
  int num_kmer = encodeKmer<int>("GGGGGGGGGGGGGG");
  EXPECT_EQ(536870911, num_kmer);
}

TEST(Kmers, IntToKmerTest) {
  std::string kmer = decodeKmer<int>(283);
  EXPECT_EQ("ACTG", kmer);
}

TEST(Kmers, LongLongToKmerTest) {
  std::string kmer = decodeKmer<long long>(119320958947904038LL);
  EXPECT_EQ("TTCGGTTCGACGTTGACCTCCATTATCT", kmer);
}

TEST(Kmers, WindowIteratorTest) {
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
