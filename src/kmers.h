#pragma once

#include <vector>
#include <string>
#include <unordered_set>

const int kNumBases = 4;
const char kBases[] = {'A', 'C', 'T', 'G'};

// Calculate move size between two kmers.
int getMove(const std::string& prev_kmer, const std::string& next_kmer);

// This function returns vector containing nextop_@dist(@kmer).
// Let Sigma = {A,C,T,G} then we define mappings
// nextop_i : Sigma^5 \rightarrow 2^{Sigma^5}.
// \forall x_1,x_2,x_3,x_4,x_5 \in Sigma, i is unsigned int:
// nextop_i(x_1 x_2 x_3 x_4 x_5) = {x_{i+1} ... x_{5} y | y \in \Sigma^i}.
std::vector<std::string> kmersInDist(const std::string& kmer, int dist);

// Returns all kmers that are in distance <= @dist.
std::unordered_set<std::string> kmersUpToDist(const std::string& kmer,
                                              int dist);

// Returns position of kmer in lexicographic order. 1-based indexing
int kmerToLexicographicPos(const std::string& kmer);

// Returns kmer that is in position @pos in lexicographic order. 1-based
// indexing.
std::string kmerInLexicographicPos(int pos, int k);

// Returns kNumBases^length. This function does not check for overflows.
inline int numKmersOf(int length);

// Converts DNA base to integer index in KBases array.
inline int baseCharToInt(char base);

// Encodes kmer to IntType. For k>14 use long long instead of int.
template <typename IntType>
IntType encodeKmer(const std::string& kmer);

template <typename IntType>
std::string decodeKmer(IntType code);

// Class that represents sliding window moving through the given string.
template <typename IntType>
class KmerWindowIterator {
 public:
  // end_window points one elements after the last element in the window.
  KmerWindowIterator(int k, const std::string::const_iterator& begin_window,
                     const std::string::const_iterator& string_end);
  bool hasNext() {
    if (end_window_ == string_end_) return false;
    return true;
  }

  IntType currentKmerCode() { return current_window_code_; }

  // Returns encoded kmer that is in the current window or -1 in case we are at
  // the end of the string.
  IntType next();

  // Return string representation of current kmer.
  std::string currentKmer() { return std::string(begin_window_, end_window_); }

 private:
  IntType most_significant_;  // kNumBases^(k-1)
  // We add one in front of the number because we want to preserve all the
  // leading zeros. Zeros represent As. See encodeKmer();
  IntType first_one_;  // kNumBases^k
  IntType current_window_code_;
  std::string::const_iterator begin_window_, end_window_, string_end_;
};

// Implementation of template classes and functions.
#include "kmers.tcc"
