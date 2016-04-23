#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "compare_samples.h"

DEFINE_string(samples_file, "",
              "First line contains ref. seq. Next line is empty and all the "
              "other lines contain samples");

DEFINE_int32(k_low, 9, "Lower bound for length of kmer.");
DEFINE_int32(step, 10, "Step size for number of samples.");
DEFINE_int32(k_upper, 30, "Upper bound for length of kmer.");

int main(int argc, char** argv) {
  google::SetUsageMessage(
      "Commandline tool for comparing kmers between Viterbi sequence and "
      "samples");
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  std::ifstream samples_file(FLAGS_samples_file);
  std::string ref;
  samples_file >> ref;

  std::string sample;
  std::vector<std::string> samples;
  while (samples_file >> sample) {
    samples.push_back(sample);
  }

  for (int k = FLAGS_k_low; k <= FLAGS_k_upper; k++) {
    for (int n_samples = 10; n_samples <= (int)samples.size();
         n_samples += 10) {
      std::random_shuffle(samples.begin(), samples.end());
      const auto& intersection_ref_size = intersectionForKmers(
          k, ref, samples.cbegin(), samples.cbegin() + n_samples);

      // Print percentage - how many kmers from ref. seq. are in the
      // intersection with samples.
      std::cout << k << "," << n_samples << ","
                << intersection_ref_size.first /
                       (double)intersection_ref_size.second << "\n";
    }
  }

  return 0;
}
