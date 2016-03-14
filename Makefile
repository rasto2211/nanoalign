CC = g++
CXXFLAGS += -g -O2 -std=gnu++11 -Wno-long-long -Wno-format -Wno-unused-result -W -Wall -Wextra -pedantic
LDFLAGS += -pthread

all: test

include google_test.mk

log2_num_test: gtest_main.a log2_num_test.o log2_num.o
hmm_test: gtest_main.a log2_num.o hmm_test.o

test: log2_num_test hmm_test

clean: 
	rm -f *.o
