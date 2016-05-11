#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <set>
#include <chrono>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "compare_samples.h"

DEFINE_string(seqs_file, "",
              "File containing ref. seq. and sequences which will be compared "
              "with ref. seq.");

DEFINE_int32(k_low, 9, "Lower bound for length of kmer.");
DEFINE_int32(k_upper, 30, "Upper bound for length of kmer.");

using std::chrono::system_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

int main(int argc, char** argv) {
  google::SetUsageMessage(
      "Commandline tool for comparing kmers between ref. sequence and other "
      "sequences. Compares individual sequences with ref.");
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  std::ifstream seqs_file(FLAGS_seqs_file);
  std::string ref;
  seqs_file >> ref;

  std::string seq;
  std::vector<std::string> seqs;
  while (seqs_file >> seq) {
    seqs.push_back(seq);
  }

  auto start = system_clock::now();
  std::cout << "k,true_positive,true_negative,false_positive,false_negative\n";
  for (int k = FLAGS_k_low; k <= FLAGS_k_upper; k++) {
    for (const StatTable& table : refVsSeqsKmers(k, ref, seqs)) {
      std::cout << k << "," << table.true_positive_ << ","
                << table.true_negative_ << "," << table.false_positive_ << ","
                << table.false_negative_ << "\n";
    }
  }
  LOG(INFO) << FLAGS_seqs_file << ": Computation of intersection of seqs took "
            << duration_cast<milliseconds>(system_clock::now() - start).count()
            << " ms";

  return 0;
}
