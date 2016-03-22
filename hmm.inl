// Implementation of templated class HMM and GaussianState.

#include <algorithm>
#include <vector>
#include <limits>
#include <stdexcept>
#include <random>
#include <chrono>

#include <cstdio>
#include <cmath>

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
    const HMM<EmissionType>::ViterbiMatrix& prob) const {
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
    Log2Num path_prob = prob[prefix][prev_state].first * transition.prob_;
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
          bestPathTo(state, prefix_len, emissions[prefix_len - 1], prob);
    }
  }

  return prob;
}

template <typename EmissionType>
std::vector<int> HMM<EmissionType>::backtrackMatrix(
    int last_state, int last_row,
    const std::function<int(int, int)>& nextState) const {
  std::vector<int> res;
  int curr_state = last_state;
  int row = last_row;
  while (row > 0) {
    res.push_back(curr_state);
    int next_state = nextState(row, curr_state);
    if (!states_[curr_state]->isSilent()) row--;
    curr_state = next_state;
  }
  res.push_back(curr_state);

  std::reverse(res.begin(), res.end());
  return res;
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

  return backtrackMatrix(
      best_terminal_state, emissions.size(),
      [&prob](int row, int state)->int { return prob[row][state].second; });
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

// Computes matrix res[i][j][k] which means:
// Sum of probabilities of all paths emiting @emissions[0...i-1] ending at
// state j and the node before j is some node which is among the first k
// predecessors in inv_transitions_[i].
template <typename EmissionType>
typename HMM<EmissionType>::ForwardMatrix HMM<EmissionType>::forwardTracking(
    const std::vector<EmissionType>& emissions) const {
  ForwardMatrix res(emissions.size() + 1,
                    std::vector<std::vector<double>>(num_states_));
  // sum_all_paths[prefix_len][state]
  // Sum of probabilities of all paths ending at @state emitting prefix of
  // emission sequence of length @prefix_len.
  std::vector<std::vector<Log2Num>> sum_all_paths(
      emissions.size() + 1, std::vector<Log2Num>(num_states_));

  // Initial values.
  for (int state = 0; state < num_states_; state++) {
    sum_all_paths[0][state] = Log2Num(0);
  }
  sum_all_paths[0][initial_state_] = Log2Num(1);

  for (int prefix_len = 1; prefix_len <= (int)emissions.size(); prefix_len++) {
    for (int state = 0; state < num_states_; state++) {
      // If the state is silent no emission is emitted. Therefore we cannot
      // extend the sequence of emission and we look at solutions with the
      // same prefix length.
      int prefix_prev;
      if (states_[state]->isSilent())
        prefix_prev = prefix_len;
      else
        prefix_prev = prefix_len - 1;

      // Sum of probabilities of all paths ending in @state and emitting
      // sequence emissions[0...prefix_prev_len-1].
      Log2Num sum = Log2Num(0);
      for (Transition transition : inv_transitions_[state]) {
        sum += transition.prob_ *
               states_[state]->prob(emissions[prefix_len - 1]) *
               sum_all_paths[prefix_prev][transition.to_state_];
        res[prefix_len][state].push_back(sum.value());
      }
      sum_all_paths[prefix_len][state] = sum;
    }
  }

  return res;
}

template <typename EmissionType>
std::vector<std::vector<int>> HMM<EmissionType>::posteriorProbSample(
    const std::vector<EmissionType>& emissions, int samples, int seed) const {
  ForwardMatrix forward_matrix = forwardTracking(emissions);

  // Sampling matrix could be calculated in forwardTracking but I split it
  // into two phases because it's better for debugging and testing. Weights in
  // discrete_distribution<int> are automatically normalized. It does not store
  // original values.
  SamplingMatrix sampling_matrix;
  for (int row = 0; row < forward_matrix.size(); row++) {
    for (int col = 0; col < forward_matrix[0].size(); col++) {
      sampling_matrix[row][col] =
          std::discrete_distribution<int>(forward_matrix[row][col]);
    }
  }

  // Weights that are used when sampling for the last state.
  std::vector<double> last_state_weights;
  const auto& last_row = forward_matrix.back();
  for (int col = 0; col < last_row.size(); col++) {
    if (last_row[col].empty()) {
      last_state_weights.push_back(0);
    } else {
      last_state_weights.push_back(last_row[col].back());
    }
  }

  std::default_random_engine generator(seed);
  std::discrete_distribution<int> last_state(last_row);
  std::vector<std::vector<int>> res(samples);
  for (int i = 0; i < samples; i++) {
    res[i] = backtrackMatrix(
        last_state(generator), sampling_matrix.size() - 1,
        [&sampling_matrix, &generator, &inv_transitions_ ](int row, int state)
                                                              ->int {
          int idx = sampling_matrix[row][state](generator);
          return inv_transitions_[state][idx].to_state_;
        });
  }
  return res;
}
