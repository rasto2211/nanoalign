#include <vector>
#include <string>
#include <stdexcept>

#include "src/pore_model.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "fast5/src/fast5.hpp"

using ::fast5::File;

TEST(PoreModel, GaussianParams1DReadScaledTest) {
  File fast5;
  fast5.open("1d_read.fast5");
  ASSERT_TRUE(fast5.is_open());
  PoreModel pore_model(fast5);

  GaussianParams params_level =
      pore_model.getGaussianParamsForLevel(kTemplate, "ACTGT");
  EXPECT_DOUBLE_EQ(1.3471538762858393, params_level.sigma_);
  EXPECT_DOUBLE_EQ(68.3454662046810597, params_level.mu_);

  GaussianParams params_sd =
      pore_model.getGaussianParamsForSd(kTemplate, "ACTGT");
  EXPECT_DOUBLE_EQ(0.1972293797204097, params_sd.sigma_);
  EXPECT_DOUBLE_EQ(0.8088667331219941, params_sd.mu_);
}

TEST(PoreModel, GaussianParams1DReadNotScaledTest) {
  File fast5;
  fast5.open("1d_read.fast5");
  ASSERT_TRUE(fast5.is_open());
  PoreModel pore_model(fast5);

  GaussianParams params_level =
      pore_model.getGaussianParamsForLevel(kTemplate, "ACTGT", false);
  EXPECT_DOUBLE_EQ(0.7669130000000000, params_level.sigma_);
  EXPECT_DOUBLE_EQ(61.8796699999999973, params_level.mu_);

  GaussianParams params_sd =
      pore_model.getGaussianParamsForSd(kTemplate, "ACTGT", false);
  EXPECT_DOUBLE_EQ(0.2803500000000000, params_sd.sigma_);
  EXPECT_DOUBLE_EQ(1.2772230000000000, params_sd.mu_);
}

TEST(PoreModel, GaussianParamsMissingStrandTest) {
  File fast5;
  fast5.open("1d_read.fast5");
  PoreModel pore_model(fast5);
  ASSERT_TRUE(fast5.is_open());
  EXPECT_THROW(pore_model.getGaussianParamsForLevel(kComplement, "ACTGT"),
               std::invalid_argument);
}

TEST(PoreModel, GaussianParamsWrongKmerTest) {
  File fast5;
  fast5.open("1d_read.fast5");
  PoreModel pore_model(fast5);
  ASSERT_TRUE(fast5.is_open());
  EXPECT_THROW(pore_model.getGaussianParamsForLevel(kTemplate, "AAAAAAA"),
               std::out_of_range);
}
