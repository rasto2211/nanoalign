#include <vector>
#include <string>
#include <set>

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
};

struct RefVsSamples {
  int samples_;
  StatTable stat_table_;
};

// Compares kmer sets of ref. seq. and every every individual seq. in input.
std::vector<StatTable> refVsSeqsKmers(int k, const std::string& ref,
                                      const std::vector<std::string>& seqs);

// Compares kmer sets of ref. seq. set of kmers for all samples.
std::vector<RefVsSamples> refVsSamplesKmers(
    int k, const std::string& ref, int step,
    const std::vector<std::string>& samples);
