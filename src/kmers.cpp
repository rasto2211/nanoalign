#include <unordered_set>

#include "kmers.h"

std::vector<std::string> kmersInDist(const std::string& kmer, int dist) {
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

int kmerToLexicographicPos(const std::string& kmer) {
  // Initialize it to number of letters lower than the last letter.
  // Number of sequences of length 1 lower than the given.
  int res = baseCharToInt(kmer.back());
  // Num of sequences of length (kmer.size()-idx-1) in the loop.
  int seqs = kNumBases;
  for (int idx = kmer.size() - 2; idx >= 0; idx--) {
    int base_idx = baseCharToInt(kmer[idx]);
    res += seqs * base_idx;
    seqs *= kNumBases;
  }

  return res + 1;
}

std::string kmerInLexicographicPos(int pos, int k) {
  int curr_pos = pos - 1;
  std::string res;
  // Number of sequences of length k-i-1.
  int seqs = numKmersOf(k - 1);
  for (int i = 0; i < k; i++) {
    res += kBases[curr_pos / seqs];
    curr_pos %= seqs;
    seqs /= kNumBases;
  }
  return res;
}

std::unordered_set<std::string> kmersUpToDist(const std::string& kmer,
                                              int dist) {
  std::unordered_set<std::string> res;
  for (int d = 0; d <= dist; d++) {
    for (const std::string& kmer_in_d : kmersInDist(kmer, d)) {
      res.insert(kmer_in_d);
    }
  }

  return res;
}
