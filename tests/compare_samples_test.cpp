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
  std::vector<long long> kmer_codes = getAllKmerCodes(3, "AACTGA");

  EXPECT_THAT(
      kmer_codes,
      ElementsAre(encodeKmer<long long>("AAC"), encodeKmer<long long>("ACT"),
                  encodeKmer<long long>("CTG"), encodeKmer<long long>("TGA")));
}

TEST(CompareSamplesTest, IntersectionForKmersTest) {
  std::pair<int, int> res = intersectionForKmers(
      3, "ACTGATAGA", {"CCTGG", "TAGTC", "GATA", "CCCTA", "CTG"});

  EXPECT_THAT(res, Pair(4, 7));
}
