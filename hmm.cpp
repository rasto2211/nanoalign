#include "hmm.h"

#include <cmath>
#include <vector>
#include <limits>

#include "log2num.h"

Log2Num GaussianState::prob(const Event& event) {
  double frac = (event.mean_ - mu_) / sigma_;
  double pi_sqrt = sqrt(2 * M_PI);
  return Log2Num((1 / (sigma_ * pi_sqrt)) * exp(-0.5 * frac * frac));
}

typedef std::pair<Log2Num, int> ProbState;

std::vector<int> HMM::runViterbi(const std::vector<Event>& events) {
  // prob[i][u] = (p,v) <=> the most probable path matching sequence
  // events[0...i-1] starting in @begin_state_ and ending in state u has
  // probability p and the state before u on this path is v.
  std::vector<std::vector<ProbState>> prob;

  // Initial probabilities.
  for (int state = 0; state < states_.size(); state++) {
    prob[0][state] = ProbState(Log2Num(0.0), kNoState);
  }
  prob[0][initial_state_] = ProbState(Log2Num(1.0), kNoState);

  for (int prefix_len = 1; prefix_len < events.size(); prefix_len++) {
    for (int state = 0; state < states_.size(); state++) {
      prob[prefix_len][state] =
          bestPathTo(state, prefix_len, events[prefix_len - 1], &prob);
    }
  }

  // Backtrack the matrix to reconstruct the best path.
  std::vector<int> res;
  int curr_state = terminal_state_;
  int prefix_len = events.size() - 1;
  while (curr_state != initial_state_) {
    res.push_back(curr_state);
    int next_state = prob[prefix_len][curr_state];
    if (!states_[curr_state].isSilent()) prefix_len--;
    curr_state = next_state;
  }
  res.push_back(curr_state);

  return res;
}

void HMM::computeInverseTransitions() {
  inv_transitions_.resize(states_.size());
  for (int state = 0; state < states_.size(); state++) {
    for (Transition transition : transitions_[state]) {
      inv_transitions_[transition.to_state_]
          .push_back({state, transition.prob_});
    }
  }
}

// Best path to @state for sequence events[0...events_prefix_len] with
// @last_event.
Log2Num HMM::bestPathTo(int state, int events_prefix_len,
                        const Event& last_event,
                        std::vector<std::vector<Log2Num>>* prob) {
  ProbState res = ProbState(Log2Num(0), kNoState);

  int prefix;
  // If the state is silent no event is emitted. Therefore the prefix for the
  // previous state is the same.
  if (states_[state].isSilent())
    prefix = events_prefix_len;
  else
    prefix = events_prefix_len - 1;

  // Try all the previous states and pick the best one.
  for (Transition transition : inv_transitions_[state]) {
    int prev_state = transition.to_state_;
    Log2Num path_prob = *prob[prefix][prev_state] * transition.prob_;
    if (res.first < path_prob) {
      res.first = path_prob;
      res.second = prev_state;
    }
  }
  res.first *= states_[state].prob(last_event) return res;
}
