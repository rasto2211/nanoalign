#include <vector>
#include <string>

#include "src/compare_samples.h"

#include "src/kmers.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using ::testing::ElementsAreArray;

TEST(CompareSamplesTest, GetHitsTest) {
  std::vector<int> num_hits =
      getNumHits(3, "ACTGTCTAG", {"ACTGACTT", "CCTGATCTCTC", "CCCCCCTAG"});
  EXPECT_THAT(num_hits, ElementsAreArray({1, 2, 0, 0, 0, 1, 1}));
}
