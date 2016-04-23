#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "compare_samples.h"

DEFINE_string(seqs_file, "",
              "File containing ref. seq. and sequences which will be compared "
              "with ref. seq.");

DEFINE_int32(k_low, 9, "Lower bound for length of kmer.");
DEFINE_int32(k_upper, 30, "Upper bound for length of kmer.");

int main(int argc, char** argv) {
  google::SetUsageMessage(
      "Commandline tool for comparing kmers between ref. sequence and other "
      "sequences. Compares individual sequences with ref.");
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  std::ifstream samples_file(FLAGS_seqs_file);
  std::string ref;
  samples_file >> ref;

  std::string seq;
  std::vector<std::string> seqs;
  while (samples_file >> seq) {
    seqs.push_back(seq);
  }

  for (const std::string& seq : seqs) {
    for (int k = FLAGS_k_low; k <= FLAGS_k_upper; k++) {
      // TODO Make convenience method for passing only one sequence to
      // intersectionForKmers(...).
      std::vector<std::string> tmp_seqs{seq};
      const auto& intersection_ref_size =
          intersectionForKmers(k, ref, tmp_seqs.cbegin(), tmp_seqs.cbegin());

      std::cout << k << "," << intersection_ref_size.first /
                                   (double)intersection_ref_size.second << "\n";
    }
  }

  return 0;
}
