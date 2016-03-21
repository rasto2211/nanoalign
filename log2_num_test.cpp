#include <iostream>

#include "log2_num.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

TEST(Log2Num, Log2NumProductTest) {
  Log2Num prod = Log2Num(0.5) * Log2Num(0.4);
  EXPECT_DOUBLE_EQ(Log2Num(0.5 * 0.4).value(), prod.value());
}

// Normal multiplication of doubles 10^-300 * 10^-300 * 10^300 would overflow.
// In this test it shouldn't.
TEST(Log2Num, Log2NumProductSmallNumberTest) {
  Log2Num prod = Log2Num(1.0e-300) * Log2Num(1.0e-300) * Log2Num(1.0e300);
  EXPECT_DOUBLE_EQ(Log2Num(1.0e-300).value(), prod.value());
}

// Test for operator *=.
TEST(Log2Num, Log2NumTimesEqualSignTest) {
  Log2Num prod = Log2Num(1.0e-300);
  prod *= Log2Num(1.0e-300);
  prod *= Log2Num(1.0e300);
  EXPECT_DOUBLE_EQ(Log2Num(1.0e-300).value(), prod.value());
}

// Test for operator +=.
TEST(Log2Num, Log2NumPlusEqualSignTest) {
  Log2Num sum = Log2Num(0.01);
  sum += Log2Num(0.1);
  sum += Log2Num(0.001);
  EXPECT_DOUBLE_EQ(Log2Num(0.111).value(), sum.value());
}

// 1.0e-300 * 0.
TEST(Log2Num, Log2NumZeroProdTest) {
  Log2Num prod = Log2Num(1.0e-300);
  prod *= Log2Num(0);
  EXPECT_DOUBLE_EQ(Log2Num(0).value(), prod.value());
}

// 1.0e-300 + 0.
TEST(Log2Num, Log2NumZeroSumTest) {
  Log2Num prod = Log2Num(1.0e-300);
  prod += Log2Num(0);
  EXPECT_DOUBLE_EQ(Log2Num(1.0e-300).value(), prod.value());
}

// Test for 0.5 + 0.1 == 0.6
TEST(Log2Num, Log2NumSumTest) {
  Log2Num sum = Log2Num(0.5) + Log2Num(0.1);
  EXPECT_DOUBLE_EQ(Log2Num(0.5 + 0.1).value(), sum.value());
}

// Test for 0 + 0.5 == 0.5
TEST(Log2Num, Log2NumSumZeroTest) {
  Log2Num sum = Log2Num(0) + Log2Num(0.5);
  EXPECT_DOUBLE_EQ(Log2Num(0.5).value(), sum.value());
}

TEST(Log2Num, Log2NumLessThanTest) { EXPECT_TRUE(Log2Num(0.5) < Log2Num(0.7)); }

TEST(Log2Num, Log2NumLessThanZerosTest) {
  EXPECT_FALSE(Log2Num(0) < Log2Num(0));
}

TEST(Log2Num, Log2NumGreaterThanOneZeroTrueTest) {
  EXPECT_TRUE(Log2Num(0.6) > Log2Num(0));
}

TEST(Log2Num, Log2NumGreaterThanZerosTest) {
  EXPECT_FALSE(Log2Num(0) > Log2Num(0));
}

TEST(Log2Num, Log2NumGreaterThanOneZeroFalseTest) {
  EXPECT_FALSE(Log2Num(0) > Log2Num(0.11));
}

TEST(Log2Num, Log2NumValue) { EXPECT_DOUBLE_EQ(0.123, Log2Num(0.123).value()); }

TEST(Log2Num, Log2NumValueZero) { EXPECT_DOUBLE_EQ(0, Log2Num(0).value()); }
