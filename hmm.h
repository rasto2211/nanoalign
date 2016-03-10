#pragma once

#include <vector>
#include <cmath>

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
  virtual Log2Num prob(const EmissionType& event) const = 0;
};

// State with no emission.
template <typename EmissionType>
class SilentState : public State<EmissionType> {
 public:
  bool isSilent() const { return true; }
  Log2Num prob(const EmissionType& /* event */) const { return Log2Num(1); }
};

// State with Gaussian emission.
class GaussianState : public State<double> {
 public:
  GaussianState(double mu, double sigma) : mu_(mu), sigma_(sigma) {}
  bool isSilent() const { return false; }
  Log2Num prob(const double& event) const;

 private:
  double mu_;
  double sigma_;
};

// Hidden Markov Model with silent states.
// It has one initial state and one terminal state.
template <typename EmissionType>
class HMM {
 public:
  HMM(int initial_state, const std::vector<State<EmissionType>*>& states,
      const std::vector<std::vector<Transition>>& transitions)
      : initial_state_(initial_state),
        num_states_(states.size()),
        states_(states),
        transitions_(transitions) {
    computeInvTransitions();
  }

  // Runs Viterbi algorithm and returns sequence of states.
  std::vector<int> runViterbi(const std::vector<EmissionType>& events) const;
  // Samples from P(state_sequence|event_sequence) and returns state sequence.
  // std::vector<int> posteriorProbSample(const std::vector<Event>& events)
  // const;

 private:
  FRIEND_TEST(HMMTest, ComputeViterbiMatrixTest);

  typedef typename std::pair<Log2Num, int> ProbStateId;
  typedef typename std::vector<std::vector<ProbStateId>> ViterbiMatrix;

  // Finds best path to @state after @steps using @prob[steps][state].
  ProbStateId bestPathTo(int state, int events_prefix_len,
                         const EmissionType& last_event,
                         ViterbiMatrix* prob) const;
  ViterbiMatrix computeViterbiMatrix(const std::vector<EmissionType>& events)
      const;
  // Computes inverse transition.
  void computeInvTransitions();

  // This constant is used in Viterbi algorithm to denote that we cannot get
  // into this state. No previous state.
  const int kNoState = -1;

  int initial_state_;
  // Number of states including the initial state.
  int num_states_;

  // List of states with emissions.
  std::vector<State<EmissionType>*> states_;
  // List of transitions from one state to another with probabilities.
  // Ids of states are from 0 to transitions_.size()-1
  std::vector<std::vector<Transition>> transitions_;
  // Inverse transitions.
  std::vector<std::vector<Transition>> inv_transitions_;
};

// Implementation of templated class.
#include "hmm.inl"
