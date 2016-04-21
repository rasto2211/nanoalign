#include <vector>
#include <string>
#include <map>

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
    int k, const std::string& ref, const std::vector<std::string>& samples);

// Returns codes of all kmers in the sequence.
std::vector<long long> getAllKmerCodes(int k, const std::string& seq);
