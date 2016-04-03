#include <vector>
#include <stdexcept>
#include <iostream>

#include <cstddef>
#include <cassert>

#include "fast5/src/fast5.hpp"

#include "pore_model.h"
#include "nano_read.h"

using ::fast5::Model_Entry;
using ::fast5::Model_Parameters;

PoreModel::PoreModel(const ::fast5::File& fast5_file) {
  for (int strand = 0; strand < 2; strand++) {
    if (!fast5_file.have_model(strand)) {
      has_strand_[strand] = false;
      continue;
    }

    has_strand_[strand] = true;
    std::vector<Model_Entry> entries = fast5_file.get_model(strand);
    params_[strand] = fast5_file.get_model_parameters(strand);

    for (Model_Entry& entry : entries) {
      const std::string& kmer(entry.kmer);
      kmer_to_model_[strand][kmer] = entry;
    }
  }
}

GaussianParams PoreModel::getGaussianParamsForLevel(int strand,
                                                    const std::string& kmer,
                                                    bool scale) const {
  if (!has_strand_[strand])
    throw std::invalid_argument("Doesn't have this strand: " +
                                std::to_string(strand));

  const Model_Entry& entry = kmer_to_model_[strand].at(kmer);
  double mu = entry.level_mean;
  double sigma = entry.level_stdv;

  if (scale) {
    mu = mu * params_[strand].scale + params_[strand].shift;
    sigma *= params_[strand].var;
  }

  return {mu, sigma};
}

GaussianParams PoreModel::getGaussianParamsForSd(int strand,
                                                 const std::string& kmer,
                                                 bool scale) const {
  if (!has_strand_[strand])
    throw std::invalid_argument("Doesn't have this strand: " +
                                std::to_string(strand));

  const Model_Entry& entry = kmer_to_model_[strand].at(kmer);
  double mu = entry.sd_mean;
  double sigma = entry.sd_stdv;

  if (scale) {
    mu *= params_[strand].scale_sd;
    sigma *= params_[strand].var_sd;
  }

  return {mu, sigma};
}
