#include <vector>
#include <string>
#include <iostream>

#include "src/hmm.h"
#include "src/move_hmm.h"
#include "src/kmers.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using ::testing::Return;
using ::testing::ElementsAreArray;
using ::testing::UnorderedElementsAreArray;
using ::testing::Contains;
using ::testing::Eq;
using ::testing::DoubleEq;
using ::testing::ContainerEq;

TEST(MoveHMMTest, ConstructEmissionsTest) {
  std::vector<GaussianParamsKmer> gaussians = {
      {"G", 1, 0.1}, {"A", 0, 0.5}, {"T", 0.5, 0.2}, {"C", 0.5, 0.1}};

  std::vector<std::unique_ptr<State<double>>> emissions =
      constructEmissions(1, gaussians);

  ASSERT_EQ(5, emissions.size());

  ASSERT_NE(nullptr, emissions[0]);
  EXPECT_EQ(SilentState<double>(), *emissions[0]);

  ASSERT_NE(nullptr, emissions[1]);
  EXPECT_EQ(GaussianState(0, 0.5), *emissions[1]);

  ASSERT_NE(nullptr, emissions[2]);
  EXPECT_EQ(GaussianState(0.5, 0.1), *emissions[2]);

  ASSERT_NE(nullptr, emissions[3]);
  EXPECT_EQ(GaussianState(0.5, 0.2), *emissions[3]);

  ASSERT_NE(nullptr, emissions[4]);
  EXPECT_EQ(GaussianState(1, 0.1), *emissions[4]);
}

const int kMoveThreshold = 3;

// total_transitions_from_state are without pseudocount.
// Pseudocount is added for every transition going from the state so we don't
// have zero probabilities.
void test_for_transition(const std::vector<std::vector<Transition>>& res,
                         const std::string& from, const std::string& to,
                         int count_with_pseudocount,
                         int total_transitions_from_state) {
  EXPECT_THAT(res[kmerToLexicographicPos(from)],
              Contains<Transition>(
                  {kmerToLexicographicPos(to),
                   Log2Num(count_with_pseudocount /
                           (double)(kmersUpToDist(from, kMoveThreshold).size() +
                                    total_transitions_from_state))}))
      << "From: " << from << "(id = " << kmerToLexicographicPos(from) << ")";
}

// Larger test with k=4.
TEST(MoveHMMTest, ConstructTransitionsLargeTest) {
  const int k = 4;
  const int kPseudoCount = 1;

  std::vector<MoveKmer> read1 = {{0, "ACTC"},
                                 {0, "ACTC"},
                                 {1, "CTCA"},
                                 {2, "CAGC"},
                                 {0, "CAGC"},
                                 {3, "CTCA"}};
  std::vector<MoveKmer> read2 = {
      {0, "CTCA"}, {1, "CAGC"}, {3, "CTCA"}, {0, "CTCA"}};

  TransitionConstructor transition_constructor(kMoveThreshold);

  int actc = kmerToLexicographicPos("ACTC");
  int ctca = kmerToLexicographicPos("CTCA");
  int cagc = kmerToLexicographicPos("CAGC");
  transition_constructor.addRead(read1);
  std::map<std::pair<int, int>, int> read1_transitions = {{{actc, actc}, 1},
                                                          {{actc, ctca}, 1},
                                                          {{ctca, cagc}, 1},
                                                          {{cagc, cagc}, 1},
                                                          {{cagc, ctca}, 1}};
  EXPECT_THAT(transition_constructor.count_for_transition_,
              ContainerEq(read1_transitions));

  transition_constructor.addRead(read2);
  std::map<std::pair<int, int>, int> read2_transitions = {{{actc, actc}, 1},
                                                          {{actc, ctca}, 1},
                                                          {{ctca, cagc}, 2},
                                                          {{ctca, ctca}, 1},
                                                          {{cagc, cagc}, 1},
                                                          {{cagc, ctca}, 2}};
  EXPECT_THAT(transition_constructor.count_for_transition_,
              ContainerEq(read2_transitions));

  std::vector<std::vector<Transition>> res =
      transition_constructor.calculateTransitions(kPseudoCount, k);

  const int kmers = 256;  // Number of kmers of length 4.
  ASSERT_EQ(kmers + 1, res.size());

  // Initial state.
  ASSERT_EQ(kmers, res[0].size());
  int id = 1;
  for (const Transition& transition : res[0]) {
    // Probability of every kmer is equal.
    EXPECT_DOUBLE_EQ(1 / (double)kmers, transition.prob_.value());
    EXPECT_EQ(id++, transition.to_state_);
  }

  // Check transitions for other states. Do not check transitions with
  // zeros count. It would be too much.
  test_for_transition(res, "ACTC", "ACTC", 2, 2);
  test_for_transition(res, "ACTC", "CTCA", 2, 2);
  test_for_transition(res, "CTCA", "CAGC", 3, 3);
  test_for_transition(res, "CTCA", "CTCA", 2, 3);
  test_for_transition(res, "CAGC", "CTCA", 3, 3);
  test_for_transition(res, "CAGC", "CAGC", 2, 3);

  // Number of transitions for every state > 0.
  for (int state = 1; state <= kmers; state++) {
    EXPECT_EQ(
        kmersUpToDist(kmerInLexicographicPos(state, k), kMoveThreshold).size(),
        res[state].size());
  }
}

TEST(MoveHMMTest, ConstructTransitionsSmallTest) {
  const int kPseudoCount = 1;
  const int k = 2;

  std::vector<MoveKmer> read1 = {
      {0, "AG"}, {1, "GA"}, {1, "AG"}, {1, "GA"}, {1, "AG"}, {2, "TG"}};

  TransitionConstructor transition_constructor(kMoveThreshold);
  transition_constructor.addRead(read1);

  int ag = kmerToLexicographicPos("AG");
  int ga = kmerToLexicographicPos("GA");
  int tg = kmerToLexicographicPos("TG");
  std::map<std::pair<int, int>, int> read1_transitions = {
      {{ag, ga}, 2}, {{ag, tg}, 1}, {{ga, ag}, 2}};
  EXPECT_THAT(transition_constructor.count_for_transition_,
              ContainerEq(read1_transitions));

  std::vector<std::vector<Transition>> res =
      transition_constructor.calculateTransitions(kPseudoCount, k);

  const int kmers = 16;  // Number of kmers of length 2.
  ASSERT_EQ(kmers + 1, res.size());

  // Initial state.
  ASSERT_EQ(kmers, res[0].size());
  int id = 1;
  for (const Transition& transition : res[0]) {
    // Probability of every kmer is equal.
    EXPECT_DOUBLE_EQ(1 / (double)kmers, transition.prob_.value());
    EXPECT_EQ(id++, transition.to_state_);
  }

  test_for_transition(res, "AG", "GA", 3, 3);
  test_for_transition(res, "AG", "TG", 2, 3);
  test_for_transition(res, "GA", "AG", 3, 2);

  // Transitions with zero count. Only pseudocount is added to achieve non-zero
  // probability. Not all zero-count transitions are listed. There are 16^2
  // transitions.
  test_for_transition(res, "AG", "AG", 1, 3);
  test_for_transition(res, "AG", "AA", 1, 3);
  test_for_transition(res, "AG", "AC", 1, 3);
  test_for_transition(res, "AG", "CC", 1, 3);
  test_for_transition(res, "AG", "TT", 1, 3);

  test_for_transition(res, "AA", "AA", 1, 0);
  test_for_transition(res, "CC", "CC", 1, 0);
  test_for_transition(res, "TT", "TT", 1, 0);

  // Number of transitions for every state > 0.
  for (int state = 1; state <= kmers; state++) {
    EXPECT_EQ(kmers, res[state].size());
  }
}

TEST(MoveHMMTest, ConstructTransitionsExceptionTest) {
  std::vector<MoveKmer> read1 = {{0, "ACG"}, {2, "GTG"}};

  const int kMoveThresholdException = 1;
  TransitionConstructor transition_constructor(kMoveThresholdException);
  EXPECT_THROW(transition_constructor.addRead(read1), std::runtime_error);
}

TEST(MoveHMMTest, StateSeqToBasesTest) {
  std::vector<std::string> kmers = {"CGTTC", "GTTCG", "TCGGA", "CGGAA",
                                    "GGAAG", "GGAAG", "GAAGT", "GAAGT",
                                    "AAGTA", "AGTAT"};
  // Convert kmers into state sequence.
  std::vector<int> states;
  for (const std::string kmer : kmers) {
    states.push_back(kmerToLexicographicPos(kmer));
  }

  EXPECT_EQ("CGTTCGGAAGTAT", stateSeqToBases(5, states));
}

TEST(MoveHMMTest, StateSeqToBasesNoStatesTest) {
  EXPECT_EQ("", stateSeqToBases(5, {}));
}
