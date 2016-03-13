#include <iostream>

#include "log2_num.h"
#include "hmm.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

TEST(GaussianStateTest, GaussianStateTest) {
  GaussianState state = GaussianState(0.5, 1.5);
  EXPECT_FALSE(state.isSilent());
  EXPECT_DOUBLE_EQ(0.2636078939238785, state.prob(0.7).value());
}

// State used only for testing purpose. Emits letters A,B,C.
class ABCState : public State<char> {
 public:
  ABCState(double a, double b, double c)
      : a_(Log2Num(a)), b_(Log2Num(b)), c_(Log2Num(c)) {}
  bool isSilent() const { return false; }
  Log2Num prob(const char& emission) const {
    if (emission == 'A')
      return a_;
    else if (emission == 'B')
      return b_;
    return c_;
  }

 private:
  Log2Num a_, b_, c_;
};

/*
** Graph of the HMM. @ means loop 2<->2.
**      +>>>>>>>>>+
**      |    @    |
** 0 -> 1 -> 2 -> 3 <-> 4
**           |          ^
**           +<<<<<>>>>>+
*/
const int kInitialState = 0;

const std::vector<State<char>*> kStates = {
    new SilentState<char>(),     new ABCState(0.3, 0.5, 0.2),
    new ABCState(0.4, 0.4, 0.2), new ABCState(0.1, 0.1, 0.8),
    new SilentState<char>()};

const std::vector<std::vector<Transition>> kTransitions = {
    {{1, Log2Num(1)}},
    {{2, Log2Num(0.7)}, {3, Log2Num(0.3)}},
    {{3, Log2Num(0.3)}, {4, Log2Num(0.3)}, {2, Log2Num(0.4)}},
    {{4, Log2Num(1)}},
    {}};

const std::vector<char> kEmissions = {'A', 'B', 'A', 'C'};

TEST(HMMTest, ComputeViterbiMatrixTest) {
  ::HMM<char> hmm = ::HMM<char>(kInitialState, kStates, kTransitions);

  std::vector<std::vector<std::pair<double, int>>> expected_matrix = {
      {{1, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}},
      {{0, -1}, {0.3, 0}, {0, -1}, {0, -1}, {0, -1}},
      {{0, -1}, {0, -1}, {0.084, 1}, {0.009, 1}, {0.0252, 2}},
      {{0, -1}, {0, -1}, {0.01344, 2}, {0.00252, 2}, {0.004032, 2}},
      {{0, -1}, {0, -1}, {0.00107520, 2}, {0.0032256, 2}, {0.0032256, 3}}};

  ::HMM<char>::ViterbiMatrix res_matrix = hmm.computeViterbiMatrix(kEmissions);

  // Check dimensions of result matrix.
  EXPECT_EQ(expected_matrix.size(), res_matrix.size());
  EXPECT_EQ(expected_matrix[0].size(), res_matrix[0].size());

  // Check that matrices are equal.
  for (int i = 0; i < (int)expected_matrix.size(); ++i) {
    for (int j = 0; j < (int)expected_matrix[0].size(); ++j) {
      EXPECT_NEAR(expected_matrix[i][j].first,
                  res_matrix[i][j].first.value(), 1.0e-15)
          << "Probability in matrix differs at (" << i << ", " << j << ").";
      EXPECT_EQ(expected_matrix[i][j].second, res_matrix[i][j].second)
          << "Previous state in matrix differs at (" << i << ", " << j << ").";
    }
  }
}
