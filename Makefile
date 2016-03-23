CC = g++
CXXFLAGS += -g -O2 -std=gnu++11 -Wno-long-long -Wno-format -Wno-unused-result -W -Wall -Wextra -pedantic
LDLIBS += -pthread

# JsonCpp
CXXFLAGS += $(shell pkg-config --cflags jsoncpp)
LDLIBS += $(shell pkg-config --libs jsoncpp)

all: tests

include google_test.mk

tests: log2_num_test hmm_test

log2_num_test: gtest_main.a log2_num_test.o log2_num.o
hmm_test: gtest_main.a log2_num.o hmm_test.o

clean: 
	rm -f *.o
	rm -f *.a
