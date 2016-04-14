#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "fast5/src/fast5.hpp"

#include "src/move_hmm.h"

DEFINE_string(list_file, "reads.txt",
              "Text file containing path to files that are going to be used "
              "for training.");

DEFINE_bool(template_strand, true,
            "Use template(true) or complement(false) strand for training.");

DEFINE_int32(move_threshold, 2,
             "Move threshold that will be used for training. It's expected "
             "that input reads do not contains greater moves otherwise "
             "exception is thrown.");

DEFINE_int32(pseudocount, 1, "Pseudocount that will be used for training.");

DEFINE_string(suffix_filename, "",
              "Suffix of the output filename. Resulting filename will be "
              "trained_move_hmm_FLAGS_suffix_filename.json");

// Kmer size.
const int k = 5;
const int kInitialState = 0;

using fast5::File;
using fast5::Event_Entry;

int main(int argc, char** argv) {
  google::SetUsageMessage("Commandline tool for training MoveHMM.");
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  Strand strand = FLAGS_template_strand ? kTemplate : kComplement;

  // Read paths to the reads from the text file.
  std::string file_path;
  std::ifstream path_list(FLAGS_list_file);
  CHECK(path_list.is_open());

  TransitionConstructor transition_constructor(FLAGS_move_threshold);
  while (path_list >> file_path) {
    try {
      File file(file_path);
      LOG(INFO) << "Processing read: " << file_path;

      if (!file.have_events(strand)) {
        LOG(ERROR) << "File " << file_path << "does not have " << strand << ".";
        continue;
      }

      // Transform Event_Entry to MoveKmer and pass it into
      // constructTransitions.
      std::vector<MoveKmer> move_kmer;
      std::vector<Event_Entry> events = file.get_events(strand);
      for (const Event_Entry& event : events) {
        move_kmer.push_back({(int)event.move, event.model_state});
      }

      transition_constructor.addRead(move_kmer);
    }
    catch (std::exception& e) {
      LOG(ERROR) << e.what();
    }
  }

  ::HMM<double> move_hmm(
      kInitialState,
      transition_constructor.calculateTransitions(FLAGS_pseudocount, k));

  std::ofstream out("trained_move_hmm_" + FLAGS_suffix_filename + ".json");
  out << move_hmm.toJsonStr();

  return 0;
}
