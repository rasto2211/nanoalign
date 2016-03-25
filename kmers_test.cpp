#include <vector>
#include <string>

#include "kmers.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using testing::StrEq;

TEST(Kmers, AllNextKmersDist0) {
  std::vector<std::string> res = allNextKmers("ACTGC", 0);
  EXPECT_THAT(res, ::testing::ElementsAreArray({StrEq("ACTGC")}));
}

TEST(Kmers, AllNextKmersDist1) {
  std::vector<std::string> res = allNextKmers("ACTGC", 1);
  EXPECT_THAT(res,
              ::testing::ElementsAreArray({StrEq("CTGCA"), StrEq("CTGCC"),
                                           StrEq("CTGCT"), StrEq("CTGCG")}));
}

TEST(Kmers, AllNextKmersDist2) {
  std::vector<std::string> res = allNextKmers("ACTGC", 2);
  EXPECT_THAT(
      res,
      ::testing::ElementsAreArray(
          {StrEq("TGCAA"), StrEq("TGCAC"), StrEq("TGCAT"), StrEq("TGCAG"),
           StrEq("TGCCA"), StrEq("TGCCC"), StrEq("TGCCT"), StrEq("TGCCG"),
           StrEq("TGCTA"), StrEq("TGCTC"), StrEq("TGCTT"), StrEq("TGCTG"),
           StrEq("TGCGA"), StrEq("TGCGC"), StrEq("TGCGT"), StrEq("TGCGG")}));
}

TEST(Kmers, AllNextKmersDistTooMuch) {
  std::vector<std::string> res = allNextKmers("AC", 4);
  EXPECT_THAT(res, ::testing::ElementsAreArray(
                       {StrEq("AA"), StrEq("AC"), StrEq("AT"), StrEq("AG"),
                        StrEq("CA"), StrEq("CC"), StrEq("CT"), StrEq("CG"),
                        StrEq("TA"), StrEq("TC"), StrEq("TT"), StrEq("TG"),
                        StrEq("GA"), StrEq("GC"), StrEq("GT"), StrEq("GG")}));
}
