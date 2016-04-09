#pragma once

#include <string>

#include <cstddef>

#include "fast5/src/fast5.hpp"

struct MoveKmer {
  int move_;
  std::string kmer_;
};

struct GaussianParamsKmer {
  std::string kmer_;
  double mu_;
  double sigma_;
};

class Fast5Reads {
 public:
  Fast5Reads(const std::string& path) : path_(path) {}
  virtual std::vector<MoveKmer> nextRead() { return std::vector<MoveKmer>(); }
  virtual bool hasNextRead() const { return false; }

 private:
  std::string path_;
};
