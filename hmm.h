#pragma once

#include <vector>
#include <cmath>

#include "log2_num.h"

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
  virtual bool isSilent() = 0;
  virtual Log2Num prob(const Event& event) = 0;
};

// State with no emission.
class SilentState : public State {
 public:
  bool isSilent() { return true; }
  Log2Num prob(const Event& event) { return Log2Num(1); }
};

// State with Gaussian emission.
class GaussianState : public State {
 public:
  bool isSilent() { return false; }
  Log2Num prob(const Event& event);
};

// Hidden Markov Model with silent states.
// It has one initial state and one terminal state.
class HMM {
 public:
  // Runs Viterbi algorithm and returns sequence of states.
  std::vector<int> runViterbi(const std::vector<Event>& events);
  // Samples from P(state_sequence|event_sequence) and returns state sequence.
  std::vector<int> posteriorProbSample(const std::vector<Event>& events);

 private:
  // Finds best path to @state after @steps using @prob[steps][state].
  Log2Num bestPathTo(int state, int events_prefix_len, const Event& last_event,
                     std::vector<std::vector<Log2Num>>* prob);
  // Computes inverse transition.
  void computeInvTransitions();

  int initial_state_;
  int terminal_state_;

  const int kNoState = -1;

  // List of states with emissions.
  std::vector<State> states_;
  // List of transitions from one state to another with probabilities.
  // Ids of states are from 0 to transitions_.size()-1
  std::vector<std::vector<Transition>> transitions_;
  // Inverse transitions.
  std::vector<std::vector<Transition>> inv_transitions_;
};
