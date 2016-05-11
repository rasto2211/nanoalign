#include <cassert>
#include <vector>
#include <string>
#include <memory>

#include <cstddef>

#include <glog/logging.h>

#include "fast5/src/fast5.hpp"

#include "move_hmm.h"
#include "kmers.h"
#include "log2_num.h"

const int kInitialState = 0;

std::vector<std::unique_ptr<State<double>>> constructEmissions(
    size_t k, const std::vector<GaussianParamsKmer>& kmer_gaussians) {
  size_t num_kmers = numKmersOf(k);
  assert(num_kmers == kmer_gaussians.size());
  std::vector<std::unique_ptr<State<double>>> res(num_kmers + 1);
  res[kInitialState] =
      std::unique_ptr<State<double>>(new SilentState<double>());
  for (const GaussianParamsKmer& gaussian : kmer_gaussians) {
    assert(gaussian.kmer_.size() == k);
    int state = kmerToLexicographicPos(gaussian.kmer_);
    res[state] = std::unique_ptr<State<double>>(
        new GaussianState(gaussian.mu_, gaussian.sigma_));
  }

  return res;
}

void TransitionConstructor::addRead(const std::vector<MoveKmer>& read) {
  if (read.empty()) return;
  // Ignore transition from initial state.
  int prev_state_id = kmerToLexicographicPos(read[0].kmer_);
  for (int i = 1; i < (int)read.size(); i++) {
    int next_state = kmerToLexicographicPos(read[i].kmer_);
    // Skip this transition if it exceeds the move threshold.
    if (read[i].move_ > move_threshold_ &&
        getSmallestMove(read[i - 1].kmer_, read[i].kmer_) > move_threshold_) {
      continue;
    }

    count_for_transition_[std::pair<int, int>(prev_state_id, next_state)]++;
    prev_state_id = next_state;
  }
}

std::vector<std::vector<Transition>>
TransitionConstructor::calculateTransitions(int pseudo_count, int k) const {
  int states = numKmersOf(k) + 1;
  // Finally, calculate transition probabilities from every state.
  std::vector<std::vector<Transition>> res(states);
  for (int id = 1; id < states; id++) {
    std::string kmer = kmerInLexicographicPos(id, k);
    std::unordered_set<std::string> next_kmers =
        kmersUpToDist(kmer, move_threshold_);

    // Number of all transitions going from @id. Sum of all counts.
    int all_transitions = 0;

    // List of transitions going from state @id.
    // Contains (next_state_id, how many times the transition occured).
    std::vector<std::pair<int, int>> transition_with_counts;
    for (const std::string& next_kmer : next_kmers) {
      int next_id = kmerToLexicographicPos(next_kmer);

      int count = pseudo_count;
      const auto& it =
          count_for_transition_.find(std::pair<int, int>(id, next_id));
      if (it != count_for_transition_.end()) {
        count += it->second;
      }

      transition_with_counts.push_back(std::pair<int, int>(next_id, count));
      all_transitions += count;
    }

    // Calculate probability for all transitions going from state @id.
    for (const std::pair<int, int>& transition : transition_with_counts) {
      double prob = transition.second / (double)all_transitions;
      res[id].push_back({transition.first, Log2Num(prob)});
    }
  }

  // Transition probabilities from initial state.
  Log2Num prob(1 / (double)(states - 1));
  for (int id = 1; id < states; id++) {
    res[kInitialState].push_back({id, prob});
  }

  return res;
}

std::string stateSeqToBases(int k, const std::vector<int>& states) {
  if (states.size() < 2) return "";

  // First state is always the initial state - silent state.
  std::string prev_kmer = kmerInLexicographicPos(states[1], k);
  std::string res(prev_kmer);
  for (int idx = 2; idx < (int)states.size(); idx++) {
    std::string next_kmer = kmerInLexicographicPos(states[idx], k);
    int move = getSmallestMove(prev_kmer, next_kmer);
    // Take suffix of length move.
    res += next_kmer.substr(next_kmer.size() - move);
    prev_kmer = next_kmer;
  }

  return res;
}
