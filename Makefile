CC = g++
WARNING_FLAGS = -Wno-long-long -Wno-format -Wno-unused-result -W -Wall -Wextra -pedantic
INC_DIRS = -Isrc/.. -Itests/.. -Ifast5/src/../..
CXXFLAGS += -g -O2 -std=gnu++11 $(WARNING_FLAGS) $(INC_DIRS)
LDLIBS += -pthread -lgflags

# JsonCpp
CXXFLAGS += $(shell pkg-config --cflags jsoncpp)
LDLIBS += $(shell pkg-config --libs jsoncpp)

# HDF5
CXXFLAGS += $(shell pkg-config --cflags hdf5)
LDLIBS += $(shell pkg-config --libs hdf5)

# libglog
CXXFLAGS += $(shell pkg-config --cflags libglog)
LDLIBS += $(shell pkg-config --libs libglog)

################################################################
#                        RULES BELOW                           #  
################################################################

all: tests tools

include tests/google_test.mk

tools: src/train_move_hmm_main src/sample_move_hmm_main src/compare_sample_kmers_main src/kmers_intersection_samples_main src/kmers_intersection_seqs_main
tests: tests/log2_num_test tests/hmm_test tests/kmers_test tests/move_hmm_test tests/compare_samples_test

src/train_move_hmm_main: src/train_move_hmm_main.o src/move_hmm.o src/kmers.o src/log2_num.o
src/sample_move_hmm_main: src/sample_move_hmm_main.o src/move_hmm.o src/kmers.o src/log2_num.o
src/compare_sample_kmers_main: src/kmers.o src/compare_samples.o
src/kmers_intersection_samples_main: src/kmers.o src/compare_samples.o
src/kmers_intersection_seqs_main: src/kmers.o src/compare_samples.o

tests/log2_num_test: tests/gtest_main.a tests/log2_num_test.o src/log2_num.o
tests/hmm_test: tests/gtest_main.a src/log2_num.o tests/hmm_test.o
tests/kmers_test: tests/gmock_main.a tests/kmers_test.o src/kmers.o
tests/pore_model_test: tests/gtest_main.a tests/pore_model_test.o src/pore_model.o
tests/move_hmm_test: tests/gmock_main.a src/move_hmm.o tests/move_hmm_test.o src/log2_num.o src/kmers.o
tests/compare_samples_test: tests/gtest_main.a src/kmers.o src/compare_samples.o

clean: 
	rm -f */*.o
