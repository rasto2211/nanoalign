#pragma once

#include <vector>
#include <cmath>
#include <memory>

#include "log2_num.h"
#include "gtest/gtest_prod.h"

// Transition from one state to another.
struct Transition {
  int to_state_;
  Log2Num prob_;
};

template <typename EmissionType>
class State {
 public:
  virtual bool isSilent() const = 0;
  virtual Log2Num prob(const EmissionType& emission) const = 0;
};

// State with no emission.
template <typename EmissionType>
class SilentState : public State<EmissionType> {
 public:
  bool isSilent() const { return true; }
  Log2Num prob(const EmissionType& /* emission */) const { return Log2Num(1); }
};

// State with Gaussian emission.
class GaussianState : public State<double> {
 public:
  GaussianState(double mu, double sigma) : mu_(mu), sigma_(sigma) {}
  bool isSilent() const { return false; }
  Log2Num prob(const double& emission) const;

 private:
  double mu_;
  double sigma_;
};

// Hidden Markov Model with silent states.
// It has one initial state.
template <typename EmissionType>
class HMM {
 public:
  // Takes ownership of State<EmissionType>*. Using unique_ptr internally.
  // States are evaluated from lowest to greatest during DP.
  // Therefore we have to put restriction on transitions going
  // to silent state. Let's say that we have transition x->y
  // and y is silent states. Then x<y has to hold true.
  // No transition can go to initial state and initial state is silent.
  HMM(int initial_state, const std::vector<State<EmissionType>*>& states,
      const std::vector<std::vector<Transition>>& transitions);

  // Runs Viterbi algorithm and returns sequence of states.
  std::vector<int> runViterbiReturnStateIds(
      const std::vector<EmissionType>& emissions) const;
  // Samples from P(state_sequence|emission_sequence) and returns state
  // sequence.
  // std::vector<int> posteriorProbSample(const std::vector<emission>&
  // emissions)
  // const;

 private:
  FRIEND_TEST(HMMTest, ComputeViterbiMatrixTest);

  typedef typename std::pair<Log2Num, int> ProbStateId;
  typedef typename std::vector<std::vector<ProbStateId>> ViterbiMatrix;
  typedef typename std::vector<std::vector<Log2Num>> Log2NumMatrix;

  // Finds best path to @state after @steps using @prob[steps][state].
  // Helper method for Viterbi algorithm.
  ProbStateId bestPathTo(int state, int emissions_prefix_len,
                         const EmissionType& last_emission,
                         const ViterbiMatrix& prob) const;
  // Computes matrix which is used in Viterbi alorithm.
  ViterbiMatrix computeViterbiMatrix(const std::vector<EmissionType>& emissions)
      const;

  // Run backward algorithm which computes sum of probabilities of all paths
  // starting in arbitrary node and emitting arbitrary suffix of emission
  // sequence.
  Log2NumMatrix backwardTracking(const std::vector<EmissionType>& emissions)
      const;
  // Helper method for backward algorithm.
  Log2Num allPathProbStartingAt(int state, int pos,
                                const EmissionType& last_emission,
                                const Log2NumMatrix& prob) const;

  // Computes inverse transition.
  void computeInvTransitions();

  // This constant is used in Viterbi algorithm to denote that we cannot get
  // into this state. No previous state.
  const int kNoState = -1;

  int initial_state_;
  // Number of states including the initial state.
  int num_states_;

  // List of states with emissions.
  std::vector<std::unique_ptr<State<EmissionType>>> states_;
  // List of transitions from one state to another with probabilities.
  // Ids of states are from 0 to transitions_.size()-1
  std::vector<std::vector<Transition>> transitions_;
  // Inverse transitions.
  std::vector<std::vector<Transition>> inv_transitions_;
};

// Implementation of templated class.
#include "hmm.inl"
