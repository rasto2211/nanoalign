#include <vector>
#include <string>
#include <queue>
#include <algorithm>

#include "kmers.h"

#define DBG(M, ...) \
  fprintf(stderr, "%s:%d: " M "\n", __FUNCTION__, __LINE__, ##__VA_ARGS__)

int baseCharToInt(char base) {
  for (int i = 0; i < kNumBases; i++) {
    if (kBases[i] == base) return i;
  }
  return -1;
}

std::vector<std::string> allNextKmers(const std::string& kmer, int dist) {
  // Cannot shift kmer by more than length of the kmer.
  dist = std::min(dist, (int)kmer.size());

  int curr_dist = 0;
  int next_dist = 1;
  std::queue<std::string> kmers_dist[2];
  kmers_dist[curr_dist].push(kmer.substr(dist));
  for (int i = 0; i < dist; i++) {
    while (!kmers_dist[curr_dist].empty()) {
      std::string curr_prefix = kmers_dist[curr_dist].front();
      kmers_dist[curr_dist].pop();
      // Append all possible bases to curr_prefix and push it into queue.
      for (int base = 0; base < kNumBases; base++) {
        kmers_dist[next_dist].push(curr_prefix + kBases[base]);
      }
    }
    std::swap(curr_dist, next_dist);
  }

  std::vector<std::string> res;
  while (!kmers_dist[curr_dist].empty()) {
    res.push_back(kmers_dist[curr_dist].front());
    kmers_dist[curr_dist].pop();
  }

  return res;
}

template <typename IntType>
IntType encodeKmer(const std::string& kmer) {
  // Numbers in base @kBases can start with AAA... which is in @kBases 000...
  // Normally, the number of leading zeros doesn't mean anything. In this case
  // we care about it. Therefore we put 1 in front of the number to not lose all
  // the zeros.
  IntType res = (IntType)1;
  for (int idx = 0; idx < (int)kmer.size(); idx++) {
    res = res * kNumBases + (IntType)baseCharToInt(kmer[idx]);
  }
  return res;
}

template <typename IntType>
std::string decodeKmer(IntType code) {
  IntType num = code;
  std::string res;
  while (num > 0) {
    res += kBases[num % kNumBases];
    num /= kNumBases;
  }
  // Get rid of 1 that we artificially added in front of the number. See
  // encodeKmer for description.
  res.pop_back();
  reverse(res.begin(), res.end());

  return res;
}

template <typename IntType>
KmerWindowIterator<IntType>::KmerWindowIterator(
    int k, const std::string::iterator& begin_window,
    const std::string::iterator& string_end)
    : begin_window_(begin_window), string_end_(string_end) {
  end_window_ = begin_window_ + k;
  current_window_code_ =
      encodeKmer<IntType>(std::string(begin_window_, end_window_));
  most_significant_ = 1;
  for (int i = 0; i < k - 1; i++) most_significant_ *= kNumBases;
  first_one_ = most_significant_ * kNumBases;
}

template <typename IntType>
IntType KmerWindowIterator<IntType>::next() {
  if (!hasNext()) return (IntType)(-1);

  IntType first_val = baseCharToInt(*begin_window_);
  IntType new_val = baseCharToInt(*end_window_);

  // Remove first char from the window.
  current_window_code_ -= most_significant_ * first_val + first_one_;
  current_window_code_ *= kNumBases;
  // Add new char to the window and add the one in front of the number.
  current_window_code_ += new_val + first_one_;

  begin_window_++;
  end_window_++;

  return current_window_code_;
}
