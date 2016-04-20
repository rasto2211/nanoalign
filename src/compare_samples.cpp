#include <vector>
#include <string>
#include <map>

#include "kmers.h"

std::vector<int> getNumHits(int k, const std::string& ref_seq,
                            const std::vector<std::string>& samples) {
  std::vector<int> ref_kmer_codes;
  KmerWindowIterator<int> kmer_window_it(k, ref_seq.begin(), ref_seq.end());
  ref_kmer_codes.push_back(kmer_window_it.currentKmerCode());
  while (kmer_window_it.hasNext()) {
    ref_kmer_codes.push_back(kmer_window_it.next());
  }

  std::vector<int> res(ref_kmer_codes.size());
  for (const std::string& sample : samples) {
    int pos = 0;
    KmerWindowIterator<int> kmer_window_it(k, sample.begin(), sample.end());
    do {
      int code = kmer_window_it.currentKmerCode();
      if (ref_kmer_codes[pos] == code) res[pos]++;
      pos++;
    } while (kmer_window_it.next() != -1 && pos < (int)res.size());
  }

  return res;
}
