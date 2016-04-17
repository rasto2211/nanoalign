// Implementation of templated class HMM and GaussianState.
#include <algorithm>
#include <vector>
#include <limits>
#include <stdexcept>
#include <random>
#include <chrono>
#include <numeric>

#include <cstdio>
#include <cmath>

#include <json/value.h>

#include "log2_num.h"

#define DBG(M, ...) \
  fprintf(stderr, "%s:%d: " M "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__)

// Parse output from macro __PRETTY_FUNCTION__ to get template type.
// The part that we are interested is "with EmissionType = char]".
template <typename T>
std::string template_to_str() {
  std::string pretty_function = __PRETTY_FUNCTION__;
  std::string emission_type = "with T = ";
  int begin = pretty_function.find(emission_type) + emission_type.size();
  int length = pretty_function.find(";", begin) - begin;
  return pretty_function.substr(begin, length);
}

inline Log2Num GaussianState::prob(const double& emission) const {
  double frac = (emission - mu_) / sigma_;
  double pi_sqrt = sqrt(2 * M_PI);
  return Log2Num((1 / (sigma_ * pi_sqrt)) * exp(-0.5 * frac * frac));
}

template <typename EmissionType>
Json::Value SilentState<EmissionType>::toJsonValue() const {
  Json::Value json_map;
  json_map["stateClass"] =
      "SilentState<" + template_to_str<EmissionType>() + ">";
  return json_map;
}

Json::Value inline GaussianState::toJsonValue() const {
  Json::Value json_map;
  json_map["stateClass"] = "GaussianState";
  json_map["params"]["mu"] = mu_;
  json_map["params"]["sigma"] = sigma_;
  return json_map;
}

template <typename EmissionType>
std::string HMM<EmissionType>::toJsonStr() const {
  Json::Value json_map;
  json_map["initial_state"] = initial_state_;

  // Serialize transitions. We need to serialize list of lists of all
  // transitions for every state.
  for (int i = 0; i < num_states_; i++) {
    const std::string& label = "transitions from " + std::to_string(i);
    json_map["transitions"][label] = Json::arrayValue;
    for (const Transition& transition : transitions_[i]) {
      Json::Value transition_json = Json::objectValue;
      transition_json["prob"] = transition.prob_.toString();
      transition_json["to_state"] = transition.to_state_;
      json_map["transitions"][label].append(transition_json);
    }
  }
  return json_map.toStyledString();
}

template <typename EmissionType>
HMM<EmissionType>::HMM(const Json::Value& hmm_json) {
  initial_state_ = hmm_json["initial_state"].asInt();
  num_states_ = hmm_json["transitions"].size();

  // Deserialize transitions.
  transitions_.resize(num_states_);
  for (int state_id = 0; state_id < num_states_; state_id++) {
    const std::string& label = "transitions from " + std::to_string(state_id);
    const Json::Value& transitions_json = hmm_json["transitions"][label];
    for (const Json::Value& trans : transitions_json) {
      transitions_[state_id].push_back(
          {trans["to_state"].asInt(), Log2Num(trans["prob"].asString())});
    }
  }
  computeInvTransitions();
}

template <typename EmissionType>
HMM<EmissionType>::HMM(int initial_state,
                       const std::vector<std::vector<Transition>>& transitions)
    : initial_state_(initial_state),
      num_states_(transitions.size()),
      // states(std::make_move_iterator(std::begin(states)),
      // std::make_move_iterator(std::end(states))),
      transitions_(transitions) {
  computeInvTransitions();
}

// These checks are done:
// 1) Initial state has to be silent.
// 2) No transitions can go to initial state.
// 3) Transition to silent state. Outgoing state has to have lower number.
template <typename EmissionType>
void HMM<EmissionType>::isValid(
    const std::vector<std::unique_ptr<State<EmissionType>>>& states) const {
  // Input validation. Checks only less expected restrictions on input.
  if (!states[initial_state_]->isSilent()) {
    throw std::invalid_argument("Initial state has to be silent.");
  }

  for (int state = 0; state < num_states_; state++) {
    for (const Transition& transition : transitions_[state]) {
      int to_state = transition.to_state_;
      if (to_state == initial_state_) {
        throw std::invalid_argument("No transitions can go to initial state.");
      }
      if (states[to_state]->isSilent() && state >= to_state) {
        throw std::invalid_argument(
            "Transition to silent state. Outgoing state has to have lower "
            "number.");
      }
    }
  }
}

// Best path to @state for sequence emissions[0...emissions_prefix_len] with
// @last_emission.
template <typename EmissionType>
typename HMM<EmissionType>::ProbStateId HMM<EmissionType>::bestPathTo(
    int state_id, const State<EmissionType>& state, int emissions_prefix_len,
    const EmissionType& last_emission,
    const HMM<EmissionType>::ViterbiMatrix& prob) const {
  ProbStateId res = ProbStateId(Log2Num(0), kNoState);

  int prefix;
  // If the state is silent no emission is emitted. Therefore the prefix for the
  // previous state is the same.
  if (state.isSilent())
    prefix = emissions_prefix_len;
  else
    prefix = emissions_prefix_len - 1;

  // Try all the previous states and pick the best one.
  for (Transition transition : inv_transitions_[state_id]) {
    int prev_state = transition.to_state_;
    Log2Num path_prob = prob[prefix][prev_state].first * transition.prob_;
    if (res.first < path_prob) {
      res.first = path_prob;
      res.second = prev_state;
    }
  }
  res.first *= state.prob(last_emission);
  return res;
}

// ViterbiMatrix[i][u] = (p,v) <=> the most probable path matching sequence
// emissions[0...i-1] starting in @begin_state_ and ending in state u has
// probability p and the state before u on this path is v.
template <typename EmissionType>
typename HMM<EmissionType>::ViterbiMatrix
HMM<EmissionType>::computeViterbiMatrix(
    const std::vector<EmissionType>& emissions,
    const std::vector<std::unique_ptr<State<EmissionType>>>& states) const {
  ViterbiMatrix prob = ViterbiMatrix(emissions.size() + 1,
                                     std::vector<ProbStateId>(num_states_));

  // Initial probabilities.
  for (int state = 0; state < num_states_; state++) {
    prob[0][state] = ProbStateId(Log2Num(0.0), kNoState);
  }
  prob[0][initial_state_] = ProbStateId(Log2Num(1.0), kNoState);

  for (int prefix_len = 1; prefix_len <= (int)emissions.size(); prefix_len++) {
    for (int state_id = 0; state_id < num_states_; state_id++) {
      prob[prefix_len][state_id] =
          bestPathTo(state_id, *states[state_id], prefix_len,
                     emissions[prefix_len - 1], prob);
    }
  }

  return prob;
}

template <typename EmissionType>
std::vector<int> HMM<EmissionType>::backtrackMatrix(
    int last_state, int last_row,
    const std::vector<std::unique_ptr<State<EmissionType>>>& states,
    const std::function<int(int, int)>& nextState) const {
  std::vector<int> res;
  int curr_state = last_state;
  int row = last_row;
  while (row > 0) {
    res.push_back(curr_state);
    int next_state = nextState(row, curr_state);
    if (!states[curr_state]->isSilent()) row--;
    curr_state = next_state;
  }
  res.push_back(curr_state);

  std::reverse(res.begin(), res.end());
  return res;
}

template <typename EmissionType>
std::vector<int> HMM<EmissionType>::runViterbiReturnStateIds(
    const std::vector<EmissionType>& emission_seq,
    const std::vector<std::unique_ptr<State<EmissionType>>>& states) const {
  // Checks is the input states and transitions are valid.
  isValid(states);

  ViterbiMatrix prob = computeViterbiMatrix(emission_seq, states);

  Log2Num best_prob = Log2Num(0);
  int best_terminal_state = 0;
  for (int i = 0; i < num_states_; ++i) {
    if (prob[emission_seq.size()][i].first > best_prob) {
      best_prob = prob[emission_seq.size()][i].first;
      best_terminal_state = i;
    }
  }

  return backtrackMatrix(
      best_terminal_state, emission_seq.size(), states,
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
// Sum of probabilities of all paths of form
// initial_state -> ... -> inv_transitions_[j][k] -> j
// and emitting emissions[0... i-1].
template <typename EmissionType>
typename HMM<EmissionType>::ForwardMatrix HMM<EmissionType>::forwardTracking(
    const std::vector<EmissionType>& emissions,
    const std::vector<std::unique_ptr<State<EmissionType>>>& states) const {
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
      if (states[state]->isSilent())
        prefix_prev = prefix_len;
      else
        prefix_prev = prefix_len - 1;

      // Sum of probabilities of all paths ending in @state and emitting
      // sequence emissions[0...prefix_prev_len-1].
      Log2Num sum = Log2Num(0);
      std::vector<Log2Num> path_probs;
      for (Transition transition : inv_transitions_[state]) {
        Log2Num path_prob = transition.prob_ *
                            states[state]->prob(emissions[prefix_len - 1]) *
                            sum_all_paths[prefix_prev][transition.to_state_];
        path_probs.push_back(path_prob);
        sum += path_prob;
      }
      sum_all_paths[prefix_len][state] = sum;

      // Normalize probabilities.
      if (!sum.isLogZero()) {
        for (const Log2Num& prob : path_probs) {
          res[prefix_len][state].push_back((prob / sum).value());
        }
      }
    }
  }

  return res;
}

template <typename EmissionType>
std::vector<std::vector<int>> HMM<EmissionType>::posteriorProbSample(
    int samples, int seed, const std::vector<EmissionType>& emission_seq,
    const std::vector<std::unique_ptr<State<EmissionType>>>& states) const {
  // Checks is the input states and transitions are valid.
  isValid(states);

  ForwardMatrix forward_matrix = forwardTracking(emission_seq, states);

  // Sampling matrix could be calculated in forwardTracking but I split it
  // into two phases because it's better for debugging and testing. Weights in
  // discrete_distribution<int> are automatically normalized. It does not store
  // original values.
  int num_states_ = states.size();
  SamplingMatrix sampling_matrix(forward_matrix.size());
  for (int row = 0; row < (int)forward_matrix.size(); row++) {
    for (int col = 0; col < num_states_; col++) {
      const auto& weights = forward_matrix[row][col];
      sampling_matrix[row].emplace_back(weights.begin(), weights.end());
    }
  }

  // Weights that are used when sampling for the last state.
  std::vector<double> last_state_weights(num_states_);
  const auto& last_row = forward_matrix.back();
  for (int col = 0; col < num_states_; col++) {
    last_state_weights[col] =
        std::accumulate(last_row[col].begin(), last_row[col].end(), 0.0);
  }

  std::default_random_engine generator(seed);
  std::discrete_distribution<int> last_state(last_state_weights.begin(),
                                             last_state_weights.end());
  std::vector<std::vector<int>> res(samples);
  for (int i = 0; i < samples; i++) {
    res[i] = backtrackMatrix(
        last_state(generator), sampling_matrix.size() - 1, states,
        [&sampling_matrix, &generator, this ](int row, int state)->int {
          int idx = sampling_matrix[row][state](generator);
          return inv_transitions_[state][idx].to_state_;
        });
  }
  return res;
}
