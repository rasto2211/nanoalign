#include <vector>
#include <string>

#include "src/compare_samples.h"

#include "src/kmers.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using ::testing::ElementsAreArray;
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
