CC = g++
CXXFLAGS += -g -O2 -std=gnu++11 -Wno-long-long -Wno-format -Wno-unused-result -W -Wall -Wextra -pedantic
LDFLAGS += -pthread

all: test

# Points to the root of Google Test, relative to where this file is.
# Remember to tweak this if you move this file.
GTEST_DIR = /usr/src/gtest

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = /usr/include/gtest/*.h \
                /usr/include/gtest/internal/*.h

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

log2_num_test: gtest_main.a log2_num_test.o log2_num.o
hmm_test: gtest_main.a log2_num.o hmm_test.o

test: log2_num_test hmm_test

clean: 
	rm -f *.o
