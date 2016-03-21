#include <iostream>
#include <stdexcept>

#include "log2_num.h"
#include "hmm.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

const double kDoubleTolerance = 1.0e-15;

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

const std::vector<std::vector<Transition>> kTransitions = {
    {{1, Log2Num(1)}},
    {{2, Log2Num(0.7)}, {3, Log2Num(0.3)}},
    {{3, Log2Num(0.3)}, {4, Log2Num(0.3)}, {2, Log2Num(0.4)}},
    {{4, Log2Num(1)}},
    {}};

const std::vector<char> kEmissions = {'A', 'B', 'A', 'C'};

std::vector<State<char>*> allocateStates() {
  return {new SilentState<char>(),     new ABCState(0.3, 0.5, 0.2),
          new ABCState(0.4, 0.4, 0.2), new ABCState(0.1, 0.1, 0.8),
          new SilentState<char>()};
}

TEST(HMMTest, ComputeInvTransitions) {
  // Inverse transitions are computed in constructor.
  ::HMM<char> hmm = ::HMM<char>(kInitialState, allocateStates(), kTransitions);

  std::vector<std::vector<Transition>> inv_transitions = hmm.inv_transitions_;

  std::vector<std::vector<Transition>> expected_inv_transitions = {
      {},
      {},
      {{1, Log2Num(0.7)}, {2, Log2Num(0.4)}},
      {{1, Log2Num(0.3)}, {2, Log2Num(0.3)}},
      {{2, Log2Num(0.3)}, {3, Log2Num(1)}}};
  // Check dimensions of result matrix.
  EXPECT_EQ(expected_inv_transitions.size(), inv_transitions.size());
  EXPECT_EQ(expected_inv_transitions[0].size(), inv_transitions[0].size());

  // Check that matrices are equal.
  for (int i = 0; i < (int)expected_inv_transitions.size(); ++i) {
    for (int j = 0; j < (int)expected_inv_transitions[i].size(); ++j) {
      EXPECT_NEAR(expected_inv_transitions[i][j].prob_.value(),
                  inv_transitions[i][j].prob_.value(), kDoubleTolerance)
          << "Probabilities differ at (" << i << ", " << j << ").";
      EXPECT_EQ(expected_inv_transitions[i][j].to_state_,
                inv_transitions[i][j].to_state_)
          << "States differ at (" << i << ", " << j << ").";
    }
  }
}

TEST(HMMTest, ComputeViterbiMatrixTest) {
  ::HMM<char> hmm = ::HMM<char>(kInitialState, allocateStates(), kTransitions);

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
                  res_matrix[i][j].first.value(), kDoubleTolerance)
          << "Probability in matrix differs at (" << i << ", " << j << ").";
      EXPECT_EQ(expected_matrix[i][j].second, res_matrix[i][j].second)
          << "Previous state in matrix differs at (" << i << ", " << j << ").";
    }
  }
}

TEST(HMMTest, RunViterbiReturnStateIdsTest) {
  HMM<char> hmm = ::HMM<char>(kInitialState, allocateStates(), kTransitions);

  std::vector<int> states = hmm.runViterbiReturnStateIds(kEmissions);
  std::vector<int> expected_states = {0, 1, 2, 2, 3};
  EXPECT_EQ(expected_states, states);
}

// When initial state is not silent constructor has to throw exception.
TEST(HMMTest, InitialStateSilentTest) {
  try {
    HMM<char> hmm =
        ::HMM<char>(kInitialState, {new ABCState(0.3, 0.5, 0.2)}, {});
    FAIL();
  }
  catch (const std::invalid_argument& err) {
    ASSERT_STREQ("Initial state has to be silent.", err.what());
  }
}

// When there's transition to initial state constructor throws exception.
TEST(HMMTest, NoTransitionToInitialStateTest) {
  try {
    HMM<char> hmm = ::HMM<char>(kInitialState, {new SilentState<char>()},
                                {{{0, Log2Num(1)}}});
    FAIL();
  }
  catch (const std::invalid_argument& err) {
    ASSERT_STREQ("No transitions can go to initial state.", err.what());
  }
}

// Let's say that we have transition x->y and y is silent states. Then x<y has
// to hold true. Otherwise exception is thrown.
TEST(HMMTest, LoopInSilentTransitionsTest) {
  try {
    std::vector<State<char>*> states = {new SilentState<char>(),
                                        new SilentState<char>(),
                                        new SilentState<char>()};
    std::vector<std::vector<Transition>> transitions = {
        {{1, Log2Num(1)}}, {{2, Log2Num(1)}}, {{1, Log2Num(1)}}};
    HMM<char> hmm = ::HMM<char>(kInitialState, states, transitions);
    FAIL();
  }
  catch (const std::invalid_argument& err) {
    ASSERT_STREQ(
        "Transition to silent state. Outgoing state has to have lower number.",
        err.what());
  }
}

TEST(HMMTest, ForwardTrackingTest) {
  ::HMM<char> hmm = ::HMM<char>(kInitialState, allocateStates(), kTransitions);

  ::HMM<char>::ForwardMatrix expected_matrix = {
      {{}, {}, {}, {}, {}},
      {{}, {0.3}, {0, 0}, {0, 0}, {0, 0}},
      {{}, {0}, {0.084, 0.084}, {0.009, 0.009}, {0.0252, 0.0342}},
      {{}, {0}, {0, 0.01344}, {0, 0.00252}, {0.004032, 0.006552}},
      {{}, {0}, {0, 0.0010752}, {0, 0.0032256}, {0.00032256, 0.00354816}}};

  HMM<char>::ForwardMatrix res_matrix = hmm.forwardTracking(kEmissions);

  for (int i = 0; i < (int)expected_matrix.size(); i++) {
    for (int j = 0; j < (int)expected_matrix[i].size(); j++) {
      for (int k = 0; k < (int)expected_matrix[i][j].size(); k++) {
        EXPECT_NEAR(res_matrix[i][j][k],
                    expected_matrix[i][j][k], kDoubleTolerance)
            << "Matrices differ at (" << i << ", " << j << ", " << k << ").";
      }
    }
  }
}
