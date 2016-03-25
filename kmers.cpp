#include <vector>
#include <string>
#include <queue>
#include <algorithm>

const int kNumBases = 4;
const char kBases[] = {'A', 'C', 'T', 'G'};

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
