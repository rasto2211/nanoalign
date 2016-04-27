#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <set>

#include "kmers.h"
#include "compare_samples.h"

#include <glog/logging.h>

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

std::set<long long> getAllKmerCodes(int k, const std::string& seq) {
  KmerWindowIterator<long long> kmer_window_it(k, seq.begin(), seq.end());

  // Return empty set in case k > seq.size().
  if (kmer_window_it.currentKmerCode() == -1) return {};

  std::set<long long> res;
  do {
    res.insert(kmer_window_it.currentKmerCode());
  } while (kmer_window_it.next() != -1);

  return res;
}

std::pair<int, int> intersectionForKmers(
    int k, const std::string& ref,
    const std::vector<std::string>::const_iterator& samples_begin,
    const std::vector<std::string>::const_iterator& samples_end) {
  std::set<long long> ref_kmers = getAllKmerCodes(k, ref);

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

long long setIntersectionSize(const std::set<long long>& s1,
                              const std::set<long long>& s2) {
  std::set<long long> smaller_set = s1;
  if (s1.size() > s2.size()) {
    smaller_set = s2;
  }

  int res = 0;
  for (long long element : smaller_set) {
    // It's set so count can be zero or one.
    res += s2.count(element);
  }

  return res;
}

StatTable calcStatsFrom(int k, int intersection_size, int ref_kmers,
                        int found_kmers) {
  long long false_positive = found_kmers - intersection_size;
  long long false_negative = ref_kmers - intersection_size;
  long long true_negative =
      numKmersOf(k) - false_negative - intersection_size - false_positive;

  return {intersection_size, true_negative, false_positive, false_negative};
}

std::vector<StatTable> refVsSeqsKmers(int k, const std::string& ref,
                                      const std::vector<std::string>& seqs) {
  std::vector<StatTable> res;
  for (const std::string& seq : seqs) {
    std::set<long long> ref_kmers = getAllKmerCodes(k, ref);
    std::set<long long> seq_kmers = getAllKmerCodes(k, seq);

    res.push_back(calcStatsFrom(k, setIntersectionSize(ref_kmers, seq_kmers),
                                ref_kmers.size(), seq_kmers.size()));
  }

  return res;
}

std::vector<StatTable> refVsSamplesKmers(
    int k, const std::string& ref, const std::vector<std::string>& samples) {
  std::vector<StatTable> res;
  std::set<long long> ref_kmers = getAllKmerCodes(k, ref);

  long long true_positive = 0;
  std::set<long long> samples_kmers_union;
  for (const std::string& sample : samples) {
    for (long long kmer_code : getAllKmerCodes(k, sample)) {
      if (samples_kmers_union.insert(kmer_code).second &&
          ref_kmers.count(kmer_code)) {
        true_positive++;  // Kmer is in the intersection.
      }
    }

    res.push_back(calcStatsFrom(k, true_positive, ref_kmers.size(),
                                samples_kmers_union.size()));
  }

  return res;
}
