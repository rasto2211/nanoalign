#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <set>

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

std::vector<long long> getAllKmerCodes(int k, const std::string& seq) {
  std::vector<long long> res;
  KmerWindowIterator<long long> kmer_window_it(k, seq.begin(), seq.end());
  do {
    res.push_back(kmer_window_it.currentKmerCode());
  } while (kmer_window_it.next() != -1);

  return res;
}

std::pair<int, int> intersectionForKmers(
    int k, const std::string& ref,
    const std::vector<std::string>::const_iterator& samples_begin,
    const std::vector<std::string>::const_iterator& samples_end) {
  std::vector<long long> ref_kmers_list = getAllKmerCodes(k, ref);
  std::set<long long> ref_kmers(ref_kmers_list.begin(), ref_kmers_list.end());

  int ref_kmers_total = ref_kmers.size();
  // Erase all kmers from ref_kmers which are in samples. All the erased kmers
  // are in intersection.
  for (auto sample_it = samples_begin; sample_it != samples_end; sample_it++) {
    for (long long kmer_code : getAllKmerCodes(k, *sample_it)) {
      ref_kmers.erase(kmer_code);
    }
  }

  return std::pair<int, int>(ref_kmers_total - ref_kmers.size(),
                             ref_kmers_total);
}
