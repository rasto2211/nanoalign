// This file contains all functions that needed to construct and work with
// MoveHMM.
#pragma once

#include <vector>
#include <string>

#include <json/value.h>

#include "hmm.h"
#include "fast5_reads.h"

// Constructs sequence of states which can be passed to HMM.
// @k - length of kmers.
// @kmer_gaussians - list of Gaussians for every kmer.
std::vector<State<double>*> constructEmissions(
    size_t k, const std::vector<GaussianParamsKmer>& kmer_gaussians);

// Constract list of transitions needed for HMM.
// @move_threshold - greatest size of move that should occur in the input. When
// this threshold is exceeded this function throws exception.
// @k - length of kmer
// @reads - All reads that are going to be used for training.
std::vector<std::vector<Transition>> constructTransitions(int move_threshold,
                                                          int pseudo_count,
                                                          int k,
                                                          Fast5Reads* reads);
