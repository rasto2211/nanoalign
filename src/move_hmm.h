// This file contains all functions that needed to construct and work with
// MoveHMM.
#pragma once

#include <vector>
#include <string>

#include <json/value.h>

#include "hmm.h"

// Represents one element in basecalled sequence.
struct MoveKmer {
  int move_;  // size of move between kmers
  std::string kmer_;
};

// Represents current level model for a given kmer.
struct GaussianParamsKmer {
  std::string kmer_;
  double mu_;
  double sigma_;
};

// DNA strand.
enum Strand {
  kTemplate = 0,
  kComplement = 1,
};

// Constructs sequence of states which can be passed to HMM.
// @k - length of kmers.
// @kmer_gaussians - list of Gaussians for every kmer.
std::vector<State<double>*> constructEmissions(
    size_t k, const std::vector<GaussianParamsKmer>& kmer_gaussians);

// This class takes reads when you call addRead() and finally constructs
// transitions when you call calculateTransitions(). Reading all reads at once
// would take too much memory so therefore it's split into two phases.
class TransitionConstructor {
 public:
  // @move_threshold - greatest size of move that should occur in the input.
  explicit TransitionConstructor(int move_threshold)
      : move_threshold_(move_threshold) {}
  // Count for every transition how many times it occurred. Results are
  // acuumulated in count_for_transition_.
  void addRead(const std::vector<MoveKmer>& read);
  // Construct list of transitions needed for HMM from count_for_transition_.
  // @k - length of kmer
  // @reads - All reads that are going to be used for training.
  std::vector<std::vector<Transition>> calculateTransitions(int pseudo_count,
                                                            int k) const;

 private:
  FRIEND_TEST(MoveHMMTest, ConstructTransitionsLargeTest);
  FRIEND_TEST(MoveHMMTest, ConstructTransitionsSmallTest);

  int move_threshold_;
  std::map<std::pair<int, int>, int> count_for_transition_;
};
