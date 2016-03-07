#include "hmm.h"

#include <cmath>
#include <vector>
#include <limits>
#include <cstdio>

#include "log2_num.h"

#define DBG(M, ...) \
  fprintf(stderr, "%s:%d: " M "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__)

Log2Num GaussianState::prob(const Event& event) const {
  double frac = (event.mean_ - mu_) / sigma_;
  double pi_sqrt = sqrt(2 * M_PI);
  return Log2Num((1 / (sigma_ * pi_sqrt)) * exp(-0.5 * frac * frac));
}

// ViterbiMatrix[i][u] = (p,v) <=> the most probable path matching sequence
// events[0...i-1] starting in @begin_state_ and ending in state u has
// probability p and the state before u on this path is v.
HMM::ViterbiMatrix HMM::computeViterbiMatrix(const std::vector<Event>& events)
    const {
  ViterbiMatrix prob;

  // Initial probabilities.
  for (int state = 0; state < num_states_; state++) {
    prob[0][state] = ProbState(Log2Num(0.0), kNoState);
  }
  prob[0][initial_state_] = ProbState(Log2Num(1.0), kNoState);

  for (int prefix_len = 1; prefix_len < (int)events.size(); prefix_len++) {
    for (int state = 0; state < num_states_; state++) {
      prob[prefix_len][state] =
          bestPathTo(state, prefix_len, events[prefix_len - 1], &prob);
    }
  }

  return prob;
}

std::vector<int> HMM::runViterbi(const std::vector<Event>& events) const {
  ViterbiMatrix prob = computeViterbiMatrix(events);

  // Backtrack the matrix to reconstruct the best path.
  std::vector<int> res;
  int curr_state = terminal_state_;
  int prefix_len = events.size() - 1;
  while (curr_state != initial_state_) {
    res.push_back(curr_state);
    int next_state = prob[prefix_len][curr_state].second;
    if (!states_[curr_state]->isSilent()) prefix_len--;
    curr_state = next_state;
  }
  res.push_back(curr_state);

  return res;
}

void HMM::computeInvTransitions() {
  inv_transitions_.resize(num_states_);
  for (int state = 0; state < num_states_; state++) {
    for (Transition transition : transitions_[state]) {
      inv_transitions_[transition.to_state_]
          .push_back({state, transition.prob_});
    }
  }
}

// Best path to @state for sequence events[0...events_prefix_len] with
// @last_event.
HMM::ProbState HMM::bestPathTo(int state, int events_prefix_len,
                               const Event& last_event,
                               HMM::ViterbiMatrix* prob) const {
  ProbState res = ProbState(Log2Num(0), kNoState);

  int prefix;
  // If the state is silent no event is emitted. Therefore the prefix for the
  // previous state is the same.
  if (states_[state]->isSilent())
    prefix = events_prefix_len;
  else
    prefix = events_prefix_len - 1;

  // Try all the previous states and pick the best one.
  for (Transition transition : inv_transitions_[state]) {
    int prev_state = transition.to_state_;
    Log2Num path_prob = (*prob)[prefix][prev_state].first * transition.prob_;
    if (res.first < path_prob) {
      res.first = path_prob;
      res.second = prev_state;
    }
  }
  res.first *= states_[state]->prob(last_event);
  return res;
}
