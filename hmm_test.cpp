#include <iostream>

#include "log2_num.h"
#include "hmm.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

TEST(GaussianStateTest, GaussianStateTest) {
  GaussianState state = GaussianState(0.5,1.5);
  EXPECT_FALSE(state.isSilent());
  EXPECT_DOUBLE_EQ(0.2636078939238785,state.prob({0.7}).value());
}

TEST(HMMTest, ComputeViterbiMatrixTest) {
  // vector<State> states = {SilentState(),GaussianState(),};
  // HMM hmm = HMM(0,5,states,transitions);
}
