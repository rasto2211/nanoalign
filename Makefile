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

################################################################
#                        RULES BELOW                           #  
################################################################

all: tests

include tests/google_test.mk

tests: tests/log2_num_test tests/hmm_test tests/kmers_test tests/pore_model_test

tests/log2_num_test: tests/gtest_main.a tests/log2_num_test.o src/log2_num.o
tests/hmm_test: tests/gtest_main.a src/log2_num.o tests/hmm_test.o
tests/kmers_test: tests/gtest_main.a tests/kmers_test.o
tests/pore_model_test: tests/gtest_main.a tests/pore_model_test.o src/pore_model.o

clean: 
	rm -f */*.o
	rm -f */*.a
