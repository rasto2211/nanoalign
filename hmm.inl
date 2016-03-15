// Implementation of templated class HMM and GaussianState.

#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>
#include <cstdio>
#include <stdexcept>

#include "log2_num.h"

#define DBG(M, ...) \
  fprintf(stderr, "%s:%d: " M "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__)

Log2Num GaussianState::prob(const double& emission) const {
  double frac = (emission - mu_) / sigma_;
  double pi_sqrt = sqrt(2 * M_PI);
  return Log2Num((1 / (sigma_ * pi_sqrt)) * exp(-0.5 * frac * frac));
}

template <typename EmissionType>
HMM<EmissionType>::HMM(int initial_state,
                       const std::vector<State<EmissionType>*>& states,
                       const std::vector<std::vector<Transition>>& transitions)
    : initial_state_(initial_state),
      num_states_(states.size()),
      states_(std::make_move_iterator(std::begin(states)),
              std::make_move_iterator(std::end(states))),
      transitions_(transitions) {
  // Input validation. Checks only less expected restrictions on input.
  if (!states_[initial_state_]->isSilent()) {
    throw std::invalid_argument("Initial state has to be silent.");
  }

  for (int state = 0; state < num_states_; state++) {
    for (const Transition& transition : transitions_[state]) {
      int to_state = transition.to_state_;
      if (to_state == initial_state_) {
        throw std::invalid_argument("No transitions can go to initial state.");
      }
      if (states_[to_state]->isSilent() && state >= to_state) {
        throw std::invalid_argument(
            "Transition to silent state. Outgoing state has to have lower "
            "number.");
      }
    }
  }

  computeInvTransitions();
}

// Best path to @state for sequence emissions[0...emissions_prefix_len] with
// @last_emission.
template <typename EmissionType>
typename HMM<EmissionType>::ProbStateId HMM<EmissionType>::bestPathTo(
    int state, int emissions_prefix_len, const EmissionType& last_emission,
    HMM<EmissionType>::ViterbiMatrix* prob) const {
  ProbStateId res = ProbStateId(Log2Num(0), kNoState);

  int prefix;
  // If the state is silent no emission is emitted. Therefore the prefix for the
  // previous state is the same.
  if (states_[state]->isSilent())
    prefix = emissions_prefix_len;
  else
    prefix = emissions_prefix_len - 1;

  // Try all the previous states and pick the best one.
  for (Transition transition : inv_transitions_[state]) {
    int prev_state = transition.to_state_;
    Log2Num path_prob = (*prob)[prefix][prev_state].first * transition.prob_;
    if (res.first < path_prob) {
      res.first = path_prob;
      res.second = prev_state;
    }
  }
  res.first *= states_[state]->prob(last_emission);
  return res;
}

// ViterbiMatrix[i][u] = (p,v) <=> the most probable path matching sequence
// emissions[0...i-1] starting in @begin_state_ and ending in state u has
// probability p and the state before u on this path is v.
template <typename EmissionType>
typename HMM<EmissionType>::ViterbiMatrix
HMM<EmissionType>::computeViterbiMatrix(
    const std::vector<EmissionType>& emissions) const {
  ViterbiMatrix prob = ViterbiMatrix(emissions.size() + 1,
                                     std::vector<ProbStateId>(num_states_));

  // Initial probabilities.
  for (int state = 0; state < num_states_; state++) {
    prob[0][state] = ProbStateId(Log2Num(0.0), kNoState);
  }
  prob[0][initial_state_] = ProbStateId(Log2Num(1.0), kNoState);

  for (int prefix_len = 1; prefix_len <= (int)emissions.size(); prefix_len++) {
    for (int state = 0; state < num_states_; state++) {
      prob[prefix_len][state] =
          bestPathTo(state, prefix_len, emissions[prefix_len - 1], &prob);
    }
  }

  return prob;
}

template <typename EmissionType>
std::vector<int> HMM<EmissionType>::runViterbiReturnStateIds(
    const std::vector<EmissionType>& emissions) const {
  ViterbiMatrix prob = computeViterbiMatrix(emissions);

  Log2Num best_prob = Log2Num(0);
  int best_terminal_state = 0;
  for (int i = 0; i < num_states_; ++i) {
    if (prob[emissions.size()][i].first > best_prob) {
      best_prob = prob[emissions.size()][i].first;
      best_terminal_state = i;
    }
  }

  // Backtrack the matrix to reconstruct the best path.
  std::vector<int> res;
  int curr_state = best_terminal_state;
  int prefix_len = emissions.size();
  while (curr_state != initial_state_) {
    res.push_back(curr_state);
    int next_state = prob[prefix_len][curr_state].second;
    if (!states_[curr_state]->isSilent()) prefix_len--;
    curr_state = next_state;
  }
  res.push_back(curr_state);
  std::reverse(res.begin(), res.end());

  return res;
}

template <typename EmissionType>
void HMM<EmissionType>::computeInvTransitions() {
  inv_transitions_.resize(num_states_);
  for (int state = 0; state < num_states_; state++) {
    for (Transition transition : transitions_[state]) {
      inv_transitions_[transition.to_state_]
          .push_back({state, transition.prob_});
    }
  }
}

template <typename EmissionType>
typename HMM<EmissionType>::Log2NumMatrix HMM<EmissionType>::backwardTracking(
    const std::vector<EmissionType>& emissions) const {
  Log2NumMatrix res(emissions.size() + 1, std::vector<Log2Num>(num_states_));

  // Initial values.
  for (int state = 0; state < num_states_; state++) {
    res[emissions.size()][state] = Log2Num(1);
  }

  for (int pos = emissions.size() - 1; pos >= 0; pos--) {
    // Attention: This loop has to be in reverse order because of silent states.
    // See description(in constructor) for the direction of edges going from the
    // silent state.
    for (int state = num_states_ - 1; state >= 0; state--) {
      // Sum of probabilities of all paths starting in @state and emitting
      // sequence emissions[pos...].
      res[pos][state] = allPathProbStartingAt(state, pos, emissions[pos], res);
    }
  }

  return res;
}

template <typename EmissionType>
Log2Num HMM<EmissionType>::allPathProbStartingAt(
    int state, int pos, const EmissionType& last_emission,
    const Log2NumMatrix& prob) const {
  Log2Num res = Log2Num(0);
  for (Transition transition : transitions_[state]) {
    int to_state = transition.to_state_;
    if (states_[state].isSilent()) {
      res += transition.prob_ * prob[pos][to_state];
    } else {
      res += transition.prob_ * states_[to_state].prob(last_emission) *
             prob[pos + 1][to_state];
    }
  }
  return res;
}
