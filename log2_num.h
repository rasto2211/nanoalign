#pragma once

#include <cmath>

// This class represents doubles in the form 2^x. Therefore only
// exponent is stored. This improves numerical stability.
// Based on:
// Mann, Tobias P. "Numerically stable hidden Markov model implementation." An
// HMM scaling tutorial (2006): 1-8.
class Log2Num {
 public:
  Log2Num() : is_log_zero_(true) {}
  // Takes number which will be converted to form 2^x.
  explicit Log2Num(double val);
  // Is it zero? Log(0) is undefined.
  bool isLogZero() const {
    return is_log_zero_;
  };
  // Sets value of the number to 2^exponent.
  void setExponent(double exponent);
  double value() const;
  Log2Num operator*(const Log2Num& num) const;
  Log2Num& operator*=(const Log2Num& num);
  Log2Num operator+(const Log2Num& num) const;
  Log2Num& operator+=(const Log2Num& num);
  bool operator<(const Log2Num& num) const;
  bool operator>(const Log2Num& num) const;

 private:
  bool is_log_zero_;
  double exponent_;
};
