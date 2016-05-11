#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <set>
#include <chrono>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "compare_samples.h"
#include "kmers.h"

DEFINE_string(samples_file, "",
              "First line contains ref. seq. Next line is empty and all the "
              "other lines contain samples");

DEFINE_int32(k_low, 9, "Lower bound for length of kmer.");
DEFINE_int32(k_upper, 30, "Upper bound for length of kmer.");

using std::chrono::system_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

int main(int argc, char** argv) {
  google::SetUsageMessage(
      "Commandline tool for comparing kmers between ref. sequence and "
      "samples. We try different sizes of random subsets of samples to compare "
      "it with reference.");
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

  auto start = system_clock::now();
  std::cout << "k,num_samples,true_positive,true_negative,false_positive,false_"
               "negative\n";
  for (int k = FLAGS_k_low; k <= FLAGS_k_upper; k++) {
    int n_samples = 0;
    for (const StatTable& stat_table : refVsSamplesKmers(k, ref, samples)) {
      n_samples++;
      std::cout << k << "," << n_samples << ", " << stat_table.true_positive_
                << "," << stat_table.true_negative_ << ","
                << stat_table.false_positive_ << ","
                << stat_table.false_negative_ << "\n";
    }
  }
  LOG(INFO) << FLAGS_samples_file
            << ": Computation of intersection of samples took "
            << duration_cast<milliseconds>(system_clock::now() - start).count()
            << " ms";

  return 0;
}
