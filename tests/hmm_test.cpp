#include <iostream>
#include <stdexcept>
#include <sstream>
#include <fstream>

#include <json/value.h>
#include <json/reader.h>

#include "src/log2_num.h"
#include "src/hmm.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

const double kDoubleTolerance = 1.0e-15;

TEST(GaussianStateTest, GaussianStateTest) {
  GaussianState state = GaussianState(0.5, 1.5);
  EXPECT_FALSE(state.isSilent());
  EXPECT_DOUBLE_EQ(0.2636078939238785, state.prob(0.7).value());
}

TEST(GaussianStateTest, GaussianStateComparisonEqTest) {
  GaussianState state1 = GaussianState(0.5, 1.5);
  GaussianState state2 = GaussianState(0.5, 1.5);
  EXPECT_TRUE(state1 == state2);
}

TEST(GaussianStateTest, GaussianStateComparisonNeqTest) {
  GaussianState state1 = GaussianState(0.5, 1.5);
  SilentState<double> state2;
  EXPECT_FALSE(state1 == state2);
}

const std::string kGaussianState =
    "{\n"
    "   \"params\" : {\n"
    "      \"mu\" : 0.123456789,\n"
    "      \"sigma\" : "
    "1\n"
    "   },\n"
    " "
    "  \"stateClass\" : \"GaussianState\"\n"
    "}\n";

TEST(GaussianStateTest, GaussianStateSerializationTest) {
  GaussianState gaussian_state(0.123456789, 1);
  EXPECT_EQ(kGaussianState, gaussian_state.toJsonStr());
}

TEST(GaussianStateTest, GaussianStateDeserializeParamsTest) {
  Json::Reader reader;
  Json::Value root;
  bool ok = reader.parse(
      "{\"mu\" : 0.123456789,\n"
      "\"sigma\" : 1\n}",
      root);
  EXPECT_TRUE(ok);
  EXPECT_EQ(0, reader.getFormattedErrorMessages().size());

  GaussianState state = GaussianState(root);
  EXPECT_DOUBLE_EQ(0.123456789, state.mu_);
  EXPECT_DOUBLE_EQ(1, state.sigma_);
}

//////////////////////////////////////////////////////////////////

TEST(SilentStateTest, SilentStateComparisonEqTest) {
  SilentState<double> state2;
  SilentState<double> state1;
  EXPECT_TRUE(state1 == state2);
}

TEST(SilentStateTest, SilentStateSerializationTest) {
  SilentState<char> silent_state;
  EXPECT_EQ(
      "{\n"
      "   \"stateClass\" : \"SilentState<char>\"\n"
      "}\n",
      silent_state.toJsonStr());
}

//////////////////////////////////////////////////////////////////

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

  Json::Value toJsonValue() const { return Json::objectValue; }

  // Fake implementation of ==. I don't need this.
  bool operator==(const State<char>& state) const {
    return typeid(*this) != typeid(state);
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

std::vector<std::unique_ptr<State<char>>> allocateStates() {
  std::vector<State<char>*> states = {
      new SilentState<char>(),     new ABCState(0.3, 0.5, 0.2),
      new ABCState(0.4, 0.4, 0.2), new ABCState(0.1, 0.1, 0.8),
      new SilentState<char>()};
  return std::vector<std::unique_ptr<State<char>>>(
      std::make_move_iterator(states.begin()),
      std::make_move_iterator(states.end()));
}

std::vector<std::unique_ptr<State<double>>> allocateDoubleStates() {
  std::vector<State<double>*> res = {
      new SilentState<double>(),   new GaussianState(0.3, 0.5),
      new GaussianState(0.4, 0.4), new GaussianState(0.1, 0.1),
      new SilentState<double>()};
  return std::vector<std::unique_ptr<State<double>>>(
      std::make_move_iterator(res.begin()), std::make_move_iterator(res.end()));
}

TEST(HMMTest, ComputeInvTransitions) {
  // Inverse transitions are computed in constructor.
  ::HMM<char> hmm = ::HMM<char>(kInitialState, kTransitions);

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
  ::HMM<char> hmm = ::HMM<char>(kInitialState, kTransitions);

  std::vector<std::vector<std::pair<double, int>>> expected_matrix = {
      {{1, -1}, {0, -1}, {0, -1}, {0, -1}, {0, -1}},
      {{0, -1}, {0.3, 0}, {0, -1}, {0, -1}, {0, -1}},
      {{0, -1}, {0, -1}, {0.084, 1}, {0.009, 1}, {0.0252, 2}},
      {{0, -1}, {0, -1}, {0.01344, 2}, {0.00252, 2}, {0.004032, 2}},
      {{0, -1}, {0, -1}, {0.00107520, 2}, {0.0032256, 2}, {0.0032256, 3}}};

  ::HMM<char>::ViterbiMatrix res_matrix =
      hmm.computeViterbiMatrix(kEmissions, allocateStates());

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
  HMM<char> hmm = ::HMM<char>(kInitialState, kTransitions);

  std::vector<int> states =
      hmm.runViterbiReturnStateIds(kEmissions, allocateStates());
  std::vector<int> expected_states = {0, 1, 2, 2, 3};
  EXPECT_EQ(expected_states, states);
}

// When initial state is not silent exception has to be thrown.
TEST(HMMTest, InitialStateSilentTest) {
  ::HMM<char> hmm = ::HMM<char>(kInitialState, {});
  std::vector<std::unique_ptr<State<char>>> states;
  states.emplace_back(new ABCState(0.3, 0.5, 0.2));

  EXPECT_THROW(hmm.runViterbiReturnStateIds(kEmissions, states),
               std::invalid_argument);

  EXPECT_THROW(hmm.posteriorProbSample(2, 0, kEmissions, states),
               std::invalid_argument);
}

// When there's transition to initial state exception has to be thrown.
TEST(HMMTest, NoTransitionToInitialStateTest) {
  ::HMM<char> hmm = ::HMM<char>(kInitialState, {{{0, Log2Num(1)}}});
  EXPECT_THROW(hmm.runViterbiReturnStateIds(kEmissions, allocateStates()),
               std::invalid_argument);
  EXPECT_THROW(hmm.posteriorProbSample(2, 0, kEmissions, allocateStates()),
               std::invalid_argument);
}

// Let's say that we have transition x->y and y is silent states. Then x<y has
// to hold true. Otherwise exception is thrown.
TEST(HMMTest, LoopInSilentTransitionsTest) {
  std::vector<std::unique_ptr<State<char>>> states;
  states.emplace_back(new SilentState<char>());
  states.emplace_back(new SilentState<char>());
  states.emplace_back(new SilentState<char>());
  std::vector<std::vector<Transition>> transitions = {
      {{1, Log2Num(1)}}, {{2, Log2Num(1)}}, {{1, Log2Num(1)}}};
  ::HMM<char> hmm = ::HMM<char>(kInitialState, transitions);
  EXPECT_THROW(hmm.runViterbiReturnStateIds(kEmissions, states),
               std::invalid_argument);
  EXPECT_THROW(hmm.posteriorProbSample(2, 0, kEmissions, states),
               std::invalid_argument);
}

TEST(HMMTest, ForwardTrackingTest) {
  ::HMM<char> hmm = ::HMM<char>(kInitialState, kTransitions);

  /*
  ** Non-normalized expected matrix:
  ** {{{}, {}, {}, {}, {}},
  ** {{}, {0.3}, {0, 0}, {0, 0}, {0, 0}},
  ** {{}, {0}, {0.084, 0}, {0.009, 0}, {0.0252, 0.009}},
  ** {{}, {0}, {0, 0.01344}, {0, 0.00252}, {0.004032, 0.00252}},
  ** {{}, {0}, {0, 0.0010752}, {0, 0.0032256}, {0.00032256, 0.0032256}}};
  */

  ::HMM<char>::ForwardMatrix expected_matrix = {
      {{}, {}, {}, {}, {}},
      {{}, {1}, {}, {}, {}},
      {{}, {}, {1, 0}, {1, 0}, {0.7368421052631578, 0.2631578947368421}},
      {{}, {}, {0, 1}, {0, 1}, {0.6153846153846153, 0.3846153846153846}},
      {{}, {}, {0, 1}, {0, 1}, {0.0909090909090909, 0.9090909090909090}}};

  HMM<char>::ForwardMatrix res_matrix =
      hmm.forwardTracking(kEmissions, allocateStates());

  EXPECT_EQ(expected_matrix.size(), res_matrix.size());
  EXPECT_EQ(expected_matrix[0].size(), res_matrix[0].size());
  for (int i = 0; i < (int)expected_matrix.size(); i++) {
    for (int j = 0; j < (int)expected_matrix[i].size(); j++) {
      EXPECT_EQ(expected_matrix[i][j].size(), res_matrix[i][j].size());
      for (int k = 0; k < (int)expected_matrix[i][j].size(); k++) {
        EXPECT_NEAR(res_matrix[i][j][k],
                    expected_matrix[i][j][k], kDoubleTolerance)
            << "Matrices differ at (" << i << ", " << j << ", " << k << ").";
      }
    }
  }
}

TEST(HMMTest, PosteriorProbSampleTest) {
  ::HMM<char> hmm = ::HMM<char>(kInitialState, kTransitions);

  std::vector<std::vector<int>> samples =
      hmm.posteriorProbSample(2, 0, kEmissions, allocateStates());

  EXPECT_THAT(samples[0][0], 0);
  EXPECT_THAT(samples[1][0], 0);
}

// Test for serialization of the whole HMM. The test json is in hmm_test.json.
TEST(HMMTest, HMMSerializationTest) {
  ::HMM<double> hmm = ::HMM<double>(kInitialState, kTransitions);
  std::ifstream json_file("hmm_test.json");
  std::stringstream ss;
  ss << json_file.rdbuf();
  EXPECT_EQ(ss.str(), hmm.toJsonStr());
}

TEST(HMMTest, HMMDeserializationTest) {
  // Convert Json file to Json::Value and then to HMM object.
  std::ifstream json_file("hmm_test.json");
  Json::CharReaderBuilder builder;
  Json::Value value;
  std::string errs;
  EXPECT_TRUE(Json::parseFromStream(builder, json_file, &value, &errs));
  ::HMM<double> hmm = ::HMM<double>(value);

  EXPECT_EQ(5, hmm.num_states_);
  EXPECT_EQ(5, hmm.transitions_.size());

  // Reads again the same file into string and serializes HMM to string and
  // compares those to strings -- original string and string from serialization.
  std::ifstream expected_json("hmm_test.json");
  std::stringstream ss;
  ss << expected_json.rdbuf();
  EXPECT_EQ(ss.str(), hmm.toJsonStr());
}
