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

// Emitted event.
struct Event {
  double mean_;
  // double duration_;
  // double std_dev_;
};

class State {
 public:
  virtual bool isSilent() const = 0;
  virtual Log2Num prob(const Event& event) const = 0;
};

// State with no emission.
class SilentState : public State {
 public:
  bool isSilent() const { return true; }
  Log2Num prob(const Event& /* event */) const { return Log2Num(1); }
};

// State with Gaussian emission.
class GaussianState : public State {
 public:
  GaussianState(double mu, double sigma) : mu_(mu), sigma_(sigma) {}
  bool isSilent() const { return false; }
  Log2Num prob(const Event& event) const;

 private:
  double mu_;
  double sigma_;
};

// Hidden Markov Model with silent states.
// It has one initial state and one terminal state.
class HMM {
 public:
  HMM(int initial_state, int terminal_state, const std::vector<State*>& states,
      const std::vector<std::vector<Transition>>& transitions)
      : initial_state_(initial_state),
        terminal_state_(terminal_state),
        states_(states),
        transitions_(transitions) {
    computeInvTransitions();
    num_states_ = states.size();
  }

  // Runs Viterbi algorithm and returns sequence of states.
  std::vector<int> runViterbi(const std::vector<Event>& events) const;
  // Samples from P(state_sequence|event_sequence) and returns state sequence.
  // std::vector<int> posteriorProbSample(const std::vector<Event>& events)
  // const;

 private:
  FRIEND_TEST(HMMTest, ComputeViterbiMatrixTest);

  typedef std::pair<Log2Num, int> ProbState;
  typedef std::vector<std::vector<ProbState>> ViterbiMatrix;

  // Finds best path to @state after @steps using @prob[steps][state].
  ProbState bestPathTo(int state, int events_prefix_len,
                       const Event& last_event, ViterbiMatrix* prob) const;
  ViterbiMatrix computeViterbiMatrix(const std::vector<Event>& events) const;
  // Computes inverse transition.
  void computeInvTransitions();

  const int kNoState = -1;

  int initial_state_;
  int terminal_state_;
  int num_states_;

  // List of states with emissions.
  std::vector<State*> states_;
  // List of transitions from one state to another with probabilities.
  // Ids of states are from 0 to transitions_.size()-1
  std::vector<std::vector<Transition>> transitions_;
  // Inverse transitions.
  std::vector<std::vector<Transition>> inv_transitions_;
};
