#include <vector>
#include <string>

#include "src/compare_samples.h"

#include "src/kmers.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using ::testing::ElementsAreArray;
using ::testing::ElementsAre;
using ::testing::Pair;

TEST(CompareSamplesTest, GetHitsAndRankTest) {
  std::vector<std::pair<int, int>> hits_rank = getNumHitsAndRank(
      3, "ACTGTCTAG", {"ACTGACTT", "CCTGATCTCTC", "CCCCCCTAG"});

  EXPECT_THAT(hits_rank,
              ElementsAreArray({Pair(1, 2), Pair(2, 0), Pair(0, 2), Pair(0, 3),
                                Pair(0, 3), Pair(1, 2), Pair(1, 1)}));
}

TEST(CompareSamplesTest, GetHitsAndRankEmptyTest) {
  std::vector<std::pair<int, int>> hits_rank = getNumHitsAndRank(3, "ACTG", {});

  EXPECT_THAT(hits_rank, ElementsAreArray({Pair(0, 0), Pair(0, 0)}));
}

TEST(CompareSamplesTest, GetAllKmersTest) {
  std::set<long long> kmer_codes = getAllKmerCodes(3, "AACTGA");

  EXPECT_THAT(
      kmer_codes,
      ElementsAre(encodeKmer<long long>("AAC"), encodeKmer<long long>("ACT"),
                  encodeKmer<long long>("CTG"), encodeKmer<long long>("TGA")));
}

TEST(CompareSamplesTest, IntersectionForKmersTest) {
  std::vector<std::string> samples = {"CCTGG", "TAGTC", "GATA", "CCCTA", "CTG"};
  std::pair<int, int> res =
      intersectionForKmers(3, "ACTGATAGA", samples.cbegin(), samples.cend());

  EXPECT_THAT(res, Pair(4, 7));
}

TEST(CompareSamplesTest, IntersectionForKmersLongLongTest) {
  std::vector<std::string> samples = {"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",
                                      "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCC"};
  std::pair<int, int> res = intersectionForKmers(
      30, "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGC", samples.cbegin(), samples.cend());

  EXPECT_THAT(res, Pair(2, 2));
}

TEST(CompareSamplesTest, RefVsSeqsKmersTest) {
  const int k = 3;
  const std::string ref = "ACTGCTA";
  const std::string seq = "AACTGCTTTT";
  std::vector<StatTable> res = refVsSeqsKmers(k, ref, {seq});
  EXPECT_EQ(1, res.size());

  StatTable stat = res[0];
  EXPECT_EQ(4, stat.true_positive_);
  EXPECT_EQ(56, stat.true_negative_);
  EXPECT_EQ(3, stat.false_positive_);
  EXPECT_EQ(1, stat.false_negative_);
}

TEST(CompareSamplesTest, RefVsSamplesKmersTest) {
  const int k = 3;
  const std::string ref = "ACTGCTA";
  std::vector<StatTable> res =
      refVsSamplesKmers(k, ref, {"ACT", "GCTTTT", "AACTGCT", "", "TTACTA"});

  // true positive, true negative, false positive, false negative.
  EXPECT_THAT(res, ElementsAreArray<StatTable>({{1, 59, 0, 4},
                                                {2, 57, 2, 3},
                                                {4, 56, 3, 1},
                                                {4, 56, 3, 1},
                                                {5, 54, 5, 0}}));
}
