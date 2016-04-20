#include <vector>
#include <string>
#include <map>

#include "kmers.h"

// Returns number of hits in samples for every kmer in ref. sequence and rank.
// Rank is number of kmers at the given position that have number of occurrences
// greater or equal than the given kmer from reference sequence.
std::vector<std::pair<int, int>> getNumHitsAndRank(
    int k, const std::string& ref_seq, const std::vector<std::string>& samples);
