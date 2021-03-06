#include <vector>
#include <string>
#include <set>
#include <iostream>

#include "kmers.h"

// Returns number of hits in samples for every kmer in ref. sequence and rank.
// Rank is number of kmers at the given position that have number of occurrences
// greater or equal than the given kmer from reference sequence.
std::vector<std::pair<int, int>> getNumHitsAndRank(
    int k, const std::string& ref_seq, const std::vector<std::string>& samples);

// Returns (intersection_size, ref_kmers_size).
// ref_kmers - set of all kmers in ref. sequence.
// samples_kmers - set of all kmers in all samples.
// intersection_size - size of intersection of these two sets.
// ref_kmers_size - |ref_kmers|
std::pair<int, int> intersectionForKmers(
    int k, const std::string& ref,
    const std::vector<std::string>::const_iterator& samples_begin,
    const std::vector<std::string>::const_iterator& samples_end);

// Returns codes of all kmers in the sequence as a sequence in order.
std::set<long long> getAllKmerCodes(int k, const std::string& seq);

struct StatTable {
  long long true_positive_;
  long long true_negative_;
  long long false_positive_;
  long long false_negative_;

  inline bool operator==(const StatTable& rhs) const {
    return (true_positive_ == rhs.true_positive_) &&
           (true_negative_ == rhs.true_negative_) &&
           (false_positive_ == rhs.false_positive_) &&
           (false_negative_ == rhs.false_negative_);
  }
};

inline std::ostream& operator<<(std::ostream& os, const StatTable& rhs) {
  os << "{" << rhs.true_positive_ << ", " << rhs.true_negative_ << ", "
     << rhs.false_positive_ << ", " << rhs.false_negative_ << "}";

  return os;
}

// Compares kmer sets of ref. seq. and every every individual seq. in input.
std::vector<StatTable> refVsSeqsKmers(int k, const std::string& ref,
                                      const std::vector<std::string>& seqs);

// Compares kmer sets of ref. seq. set of kmers for all samples.
std::vector<StatTable> refVsSamplesKmers(
    int k, const std::string& ref, const std::vector<std::string>& samples);
