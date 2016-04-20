#include <vector>
#include <string>
#include <map>
#include <algorithm>

#include "kmers.h"

std::vector<std::pair<int, int>> getNumHitsAndRank(
    int k, const std::string& ref_seq,
    const std::vector<std::string>& samples) {
  std::vector<int> ref_kmer_codes;
  KmerWindowIterator<int> kmer_window_it(k, ref_seq.begin(), ref_seq.end());
  ref_kmer_codes.push_back(kmer_window_it.currentKmerCode());
  while (kmer_window_it.hasNext()) {
    ref_kmer_codes.push_back(kmer_window_it.next());
  }

  // kmers_at_position[position][kmer_code] - number of occurrences of
  // @kmer_code at @position.
  std::vector<std::map<int, int>> kmers_at_position(ref_kmer_codes.size());
  for (const std::string& sample : samples) {
    int pos = 0;
    KmerWindowIterator<int> kmer_window_it(k, sample.begin(), sample.end());
    do {
      kmers_at_position[pos][kmer_window_it.currentKmerCode()]++;
      pos++;
    } while (kmer_window_it.next() != -1 && pos < (int)ref_kmer_codes.size());
  }

  std::vector<std::pair<int, int>> res(ref_kmer_codes.size());
  for (int pos = 0; pos < (int)res.size(); pos++) {
    int ref_kmer_code = ref_kmer_codes[pos];
    int ref_kmer_count = kmers_at_position[pos][ref_kmer_code];
    int rank = 0;
    // Rank is number of kmers which have count greater or equal than count of
    // reference kmer.
    for (const std::pair<int, int>& kmer_count : kmers_at_position[pos]) {
      if (kmer_count.first != ref_kmer_code &&
          kmer_count.second >= ref_kmer_count) {
        rank++;
      }
    }
    res[pos] = std::pair<int, int>(ref_kmer_count, rank);
  }

  return res;
}
