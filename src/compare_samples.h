#include <vector>
#include <string>
#include <map>

#include "kmers.h"

// Returns number of hits in samples for every kmer in ref. sequence.
std::vector<int> getNumHits(int k, const std::string& ref_seq,
                            const std::vector<std::string>& samples);
