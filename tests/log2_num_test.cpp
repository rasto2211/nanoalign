#include <iostream>

#include "src/log2_num.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

TEST(Log2NumTest, Log2NumProductTest) {
  Log2Num prod = Log2Num(0.5) * Log2Num(0.4);
  EXPECT_DOUBLE_EQ(Log2Num(0.5 * 0.4).value(), prod.value());
}

// Normal multiplication of doubles 10^-300 * 10^-300 * 10^300 would overflow.
// In this test it shouldn't.
TEST(Log2NumTest, Log2NumProductSmallNumberTest) {
  Log2Num prod = Log2Num(1.0e-300) * Log2Num(1.0e-300) * Log2Num(1.0e300);
  EXPECT_DOUBLE_EQ(Log2Num(1.0e-300).value(), prod.value());
}

// Test for operator *=.
TEST(Log2NumTest, Log2NumTimesEqualSignTest) {
  Log2Num prod = Log2Num(1.0e-300);
  prod *= Log2Num(1.0e-300);
  prod *= Log2Num(1.0e300);
  EXPECT_DOUBLE_EQ(Log2Num(1.0e-300).value(), prod.value());
}

// Test for operator +=.
TEST(Log2NumTest, Log2NumPlusEqualSignTest) {
  Log2Num sum = Log2Num(0.01);
  sum += Log2Num(0.1);
  sum += Log2Num(0.001);
  EXPECT_DOUBLE_EQ(Log2Num(0.111).value(), sum.value());
}

// 1.0e-300 * 0.
TEST(Log2NumTest, Log2NumZeroProdTest) {
  Log2Num prod = Log2Num(1.0e-300);
  prod *= Log2Num(0);
  EXPECT_DOUBLE_EQ(Log2Num(0).value(), prod.value());
}

// 1.0e-300 + 0.
TEST(Log2NumTest, Log2NumZeroSumTest) {
  Log2Num prod = Log2Num(1.0e-300);
  prod += Log2Num(0);
  EXPECT_DOUBLE_EQ(Log2Num(1.0e-300).value(), prod.value());
}

// Test for 0.5 + 0.1 == 0.6
TEST(Log2NumTest, Log2NumSumTest) {
  Log2Num sum = Log2Num(0.5) + Log2Num(0.1);
  EXPECT_DOUBLE_EQ(Log2Num(0.5 + 0.1).value(), sum.value());
}

// Test for 0 + 0.5 == 0.5
TEST(Log2NumTest, Log2NumSumZeroTest) {
  Log2Num sum = Log2Num(0) + Log2Num(0.5);
  EXPECT_DOUBLE_EQ(Log2Num(0.5).value(), sum.value());
}

TEST(Log2NumTest, Log2NumLessThanTest) {
  EXPECT_TRUE(Log2Num(0.5) < Log2Num(0.7));
}

TEST(Log2NumTest, Log2NumLessThanZerosTest) {
  EXPECT_FALSE(Log2Num(0) < Log2Num(0));
}

TEST(Log2NumTest, Log2NumGreaterThanOneZeroTrueTest) {
  EXPECT_TRUE(Log2Num(0.6) > Log2Num(0));
}

TEST(Log2NumTest, Log2NumGreaterThanZerosTest) {
  EXPECT_FALSE(Log2Num(0) > Log2Num(0));
}

TEST(Log2NumTest, Log2NumGreaterThanOneZeroFalseTest) {
  EXPECT_FALSE(Log2Num(0) > Log2Num(0.11));
}

TEST(Log2NumTest, Log2NumValueTest) {
  EXPECT_DOUBLE_EQ(0.123, Log2Num(0.123).value());
}

TEST(Log2NumTest, Log2NumValueZeroTest) {
  EXPECT_DOUBLE_EQ(0, Log2Num(0).value());
}

TEST(Log2NumTest, SerializationTest) {
  std::ostringstream os;
  os << Log2Num(1.234e-5);
  EXPECT_EQ("2^-16.306298079949482", os.str());
}

TEST(Log2NumTest, SerializationZeroTest) {
  std::ostringstream os;
  os << Log2Num(0);
  EXPECT_EQ("LOG2_ZERO", os.str());
}

TEST(Log2NumTest, DeserializationTest) {
  std::istringstream is("2^-16.306298079949482");
  Log2Num num;
  is >> num;
  EXPECT_DOUBLE_EQ(1.234e-5, num.value());
}

TEST(Log2NumTest, DeserializationZeroTest) {
  std::istringstream is("LOG2_ZERO");
  Log2Num num;
  is >> num;
  EXPECT_DOUBLE_EQ(0, num.value());
}

TEST(Log2NumTest, toStringTest) {
  Log2Num num(1.234e-5);
  EXPECT_EQ("2^-16.306298079949482", num.toString());
}

TEST(Log2NumTest, ConstructLog2NumFromStringTest) {
  Log2Num num("2^-16.306298079949482");
  EXPECT_DOUBLE_EQ(1.234e-5, num.value());
}

TEST(Log2NumTest, EqualityTest) {
  EXPECT_TRUE(Log2Num(1.0e-16) == Log2Num(2.0e-16));
  EXPECT_TRUE(Log2Num(1.0e-15) == Log2Num(1.0e-15 + 1.0e-16));
  EXPECT_FALSE(Log2Num(1.0e-14) == Log2Num(2.0e-14));
}

TEST(Log2NumTest, NonEqualityTest) {
  EXPECT_FALSE(Log2Num(1.0e-16) != Log2Num(2.0e-16));
  EXPECT_FALSE(Log2Num(1.0e-15) != Log2Num(1.0e-15 + 1.0e-16));
  EXPECT_TRUE(Log2Num(1.0e-14) != Log2Num(2.0e-14));
}

TEST(Log2NumTest, DivisionTest) {
  Log2Num a(0.5);
  Log2Num b(0.4);
  b /= a;
  EXPECT_EQ(0.8, b.value());
  EXPECT_EQ(0.8, (Log2Num(0.4) / Log2Num(0.5)).value());
}

TEST(Log2NumTest, DivisionZeroTest) {
  Log2Num a(0.5);
  Log2Num b(0);
  b /= a;
  EXPECT_EQ(0, b.value());
  EXPECT_EQ(0, (Log2Num(0) / Log2Num(0.5)).value());
}
