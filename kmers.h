#pragma once

#include <vector>
#include <string>

const int kNumBases = 4;
const char kBases[] = {'A', 'C', 'T', 'G'};

// This function returns vector containing nextop_@dist(@kmer).
// Let Sigma = {A,C,T,G} then we define mappings
// nextop_i : Sigma^5 \rightarrow 2^{Sigma^5}.
// \forall x_1,x_2,x_3,x_4,x_5 \in Sigma, i is unsigned int:
// nextop_i(x_1 x_2 x_3 x_4 x_5) = {x_{i+1} ... x_{5} y | y \in \Sigma^i}.
std::vector<std::string> allNextKmers(const std::string& kmer, int dist);

// Converts DNA base to integer index in KBases array.
int baseCharToInt(char base);

// Encodes kmer to IntType. For k>14 use long long instead of int.
template <typename IntType>
IntType encodeKmer(const std::string& kmer);

template <typename IntType>
std::string decodeKmer(IntType code);

/*
template <typename IntType>
IntType moveKmerWindow(int, IntType this_window_code,
                       std::string::iterator *begin_window,
                       std::string::iterator *end_window,
                       std::string::iterator *string_end) {
    if (end_window == string_end) return (IntType)-1;

    IntType first_val = **begin_window;
    IntType new_val = **end_window;
}
*/

// Class that represents sliding window moving through the given string.
template <typename IntType>
class KmerWindowIterator {
 public:
  // end_window points one elements after the last element in the window.
  KmerWindowIterator(int k, const std::string::iterator& begin_window,
                     const std::string::iterator& end_window,
                     const std::string::iterator& string_end)
      : k_(k),
        begin_window_(begin_window),
        end_window_(end_window),
        string_end_(string_end) {}

  bool hasNext() {
    if (end_window_ == string_end_) return true;
    return false;
  }
  // Returns encoded kmer that is in the current window or -1 in case we are at
  // the end of the string.
  IntType next();
  // Return string representation of current kmer.
  std::string currentKmer() { return }

 private:
  int k_;  // Length of kmer.
  std::string::iterator begin_window_, end_window_, string_end_;
};

#include "kmers.inl"
