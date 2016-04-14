
#include <cstddef>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cstddef>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "fast5/src/fast5.hpp"

#include "src/move_hmm.h"
#include "src/model_params_corrections.h"

#include <json/value.h>
#include <json/reader.h>

DEFINE_string(list_file, "reads.txt",
              "Text file containing path to files that are going to be used "
              "for training.");

DEFINE_bool(template_strand, true,
            "Use template(true) or complement(false) strand for training.");

DEFINE_string(trained_move_hmm, "",
              "Path to JSON file containing serialized MoveHMM.");

DEFINE_int32(samples, 100, "Number of samples.");

using ::fast5::File;
using ::fast5::Event_Entry;
using ::fast5::Model_Entry;
using ::fast5::Model_Parameters;

const int k = 5;  // length of kmer

std::string getFilenameFrom(const std::string& path) {
  size_t last_slash = path.find_last_of('/');
  if (last_slash == std::string::npos) return path;
  return path.substr(last_slash + 1);
}

int main(int argc, char** argv) {
  google::SetUsageMessage(
      "Commandline tool for sampling from posterior probability of MoveHMM.");
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  Strand strand = FLAGS_template_strand ? kTemplate : kComplement;

  std::string file_path;
  std::ifstream path_list(FLAGS_list_file);
  CHECK(path_list.is_open());

  // Parse json file with MoveHMM.
  Json::Value value;
  std::string errs;
  std::ifstream json_file(FLAGS_trained_move_hmm);
  Json::CharReaderBuilder builder;
  CHECK(Json::parseFromStream(builder, json_file, &value, &errs)) << errs;
  ::HMM<double> hmm = ::HMM<double>(value);

  srand(time(0));
  while (path_list >> file_path) {
    try {
      File file(file_path);
      LOG(INFO) << "Processing read: " << file_path;

      if (!file.have_events(strand)) {
        LOG(ERROR) << "File " << file_path << "does not have " << strand << ".";
        continue;
      }
      if (!file.have_model(strand)) {
        LOG(ERROR) << "File " << file_path << "does not have model for "
                   << strand << ".";
        continue;
      }

      // Get current levels.
      std::vector<double> current_levels;
      std::vector<Event_Entry> events = file.get_events(strand);
      for (const Event_Entry& event : events) {
        current_levels.push_back(event.mean);
      }

      // Construct states for given HMM.
      std::vector<Model_Entry> kmer_models = file.get_model(strand);
      Model_Parameters model_params = file.get_model_parameters(strand);
      std::vector<GaussianParamsKmer> gaussian_kmer;
      for (const Model_Entry& model_entry : kmer_models) {
        Gaussian scaled_gaussian = scaleGaussianCurrentLevel(
            {model_entry.level_mean, model_entry.level_stdv}, model_params);
        gaussian_kmer.push_back(
            {model_entry.kmer, scaled_gaussian.mu_, scaled_gaussian.sigma_});
      }
      const std::vector<State<double>*>& states =
          constructEmissions(k, gaussian_kmer);

      // Replace .fast5 with .samples extension. That'll be the output file.
      std::string filename = getFilenameFrom(file_path);
      int extension_pos = filename.find_last_of('.');
      std::string out_filename =
          filename.replace(extension_pos + 1, 5, "samples");

      std::ofstream out_file;

      // Run Viterbi algorithm.
      std::vector<int> viterbi_seq =
          hmm.runViterbiReturnStateIds(current_levels, states);
      out_file << stateSeqToBases(k, viterbi_seq) << "\n\n";

      // Sample from posterior probability.
      int seed = rand();
      std::vector<std::vector<int>> samples =
          hmm.posteriorProbSample(FLAGS_samples, seed, current_levels, states);
      for (const auto& sample : samples) {
        out_file << stateSeqToBases(k, sample) << "\n";
      }
    }
    catch (std::exception& e) {
      LOG(ERROR) << e.what();
    }
  }

  return 0;
}
