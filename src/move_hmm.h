#pragma once

#include <vector>
#include <string>

#include <json/value.h>

#include "hmm.h"
#include "fast5_reads.h"

std::vector<State<double>*> constructEmissions(
    size_t k, const std::vector<GaussianParamsKmer>& kmer_gaussians);

std::vector<std::vector<Transition>> constructTransitions(int move_threshold,
                                                          int pseudo_count,
                                                          int k,
                                                          Fast5Reads* reads);
