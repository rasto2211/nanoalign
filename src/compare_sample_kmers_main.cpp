// Commandline tool for comparing kmers between reference read and samples.

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <dirent.h>
#include <cstdio>
#include <cstring>

#include <gflags/gflags.h>
#include <glog/logging.h>

// First line of every input file contains reference sequence. Next line is
// empty and afterwards there couple lines with samples.
DEFINE_string(reads_dir, "",
              "Full path to directory containing reads with samples");

int main(int argc, char **argv) {
  google::SetUsageMessage(
      "Commandline tool for comparing kmers between Viterbi sequence and "
      "samples");
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  DIR *dir = opendir(FLAGS_reads_dir.c_str());
  if (dir == nullptr) {
    LOG(FATAL) << "Error opening directory " << FLAGS_reads_dir << ": "
               << strerror(errno);
  }

  struct dirent *dir_entry;
  while ((dir_entry = readdir(dir)) != nullptr) {
    std::string filename(dir_entry->d_name);
    std::string extension(".samples");

    if (filename.size() <= extension.size() ||
        filename.compare(filename.size() - extension.size(),
                         extension.size(), extension) != 0) {
      continue;
    }

    std::ifstream samples_file(filename);
    LOG(INFO) << "Reading file: " << filename;

    std::string ref_seq;
    samples_file >> ref_seq;
    std::string line;
    samples_file >> line;
    LOG(INFO) << line;
  }
  closedir(dir);

  return 0;
}
