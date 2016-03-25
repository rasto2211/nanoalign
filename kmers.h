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

// Converts DNA base to integer in KBases array.
int baseCharToInt(char base);
