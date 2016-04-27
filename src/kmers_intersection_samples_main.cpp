#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <set>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "compare_samples.h"
#include "kmers.h"

DEFINE_string(samples_file, "",
              "First line contains ref. seq. Next line is empty and all the "
              "other lines contain samples");

DEFINE_int32(k_low, 9, "Lower bound for length of kmer.");
DEFINE_int32(step, 10, "Step size for number of samples.");
DEFINE_int32(k_upper, 30, "Upper bound for length of kmer.");

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

  for (int k = FLAGS_k_low; k <= FLAGS_k_upper; k++) {
    for (const RefVsSamples& ref_vs_samples :
         refVsSamplesKmers(k, ref, FLAGS_step, samples)) {
      const StatTable table = ref_vs_samples.stat_table_;
      std::cout << k << "," << ref_vs_samples.samples_ << ","
                << table.true_positive_ << "," << table.true_negative_ << ","
                << table.false_positive_ << "," << table.false_negative_;
    }
  }

  return 0;
}
