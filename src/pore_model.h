#pragma once

#include <vector>
#include <unordered_map>

#include <cassert>

#include "fast5/src/fast5.hpp"

struct GaussianParams {
  double mu_;
  double sigma_;
};

enum Strand {
  kTemplate = 0,
  kComplement = 1,
};

// Model for both nanopore strands.
class PoreModel {
 public:
  PoreModel(const ::fast5::File& fast5_file);
  // Returns GaussianParams for current level - mean current.
  GaussianParams getGaussianParamsForLevel(Strand strand,
                                           const std::string& kmer,
                                           bool scale = true) const;
  // Return GaussianParams for standard deviation of current for kmer.
  GaussianParams getGaussianParamsForSd(Strand strand, const std::string& kmer,
                                        bool scale = true) const;

 private:
  // Some reads might not have both strands.
  bool has_strand_[2];
  // entries_[strand][code] stores Model_Entry for kmer with @code in @strand.
  // TODO - could be more efficient - use lexicographic order.
  std::unordered_map<std::string, ::fast5::Model_Entry> kmer_to_model_[2];
  // Scaling and other parameters for Gaussian models.
  ::fast5::Model_Parameters params_[2];
};