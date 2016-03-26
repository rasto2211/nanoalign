#pragma once

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

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
  // Takes string in form 2^exponent and constructs Log2Num from it.
  explicit Log2Num(const std::string& val_str) {
    std::istringstream is(val_str);
    is >> *this;
  }
  // Is it zero? Log(0) is undefined.
  bool isLogZero() const {
    return is_log_zero_;
  };
  // Sets value of the number to 2^exponent.
  void setExponent(double exponent);
  double value() const;
  // Log2Num is written in the form 2^exponent to string.
  std::string toString() const {
    std::ostringstream os;
    os << *this;
    return os.str();
  }

  Log2Num operator*(const Log2Num& num) const;
  Log2Num& operator*=(const Log2Num& num);
  Log2Num operator+(const Log2Num& num) const;
  Log2Num& operator+=(const Log2Num& num);
  bool operator<(const Log2Num& num) const;
  bool operator>(const Log2Num& num) const;

  friend inline std::ostream& operator<<(std::ostream& os, const Log2Num& num);
  friend inline std::istream& operator>>(std::istream& is, Log2Num& num);

 private:
  bool is_log_zero_;
  double exponent_;
};

inline std::ostream& operator<<(std::ostream& os, const Log2Num& num) {
  if (num.isLogZero()) {
    os << "LOG2_ZERO";
  } else {
    os.precision(std::numeric_limits<double>::max_digits10);
    os << "2^" << num.exponent_;
  }
  return os;
}

inline std::istream& operator>>(std::istream& is, Log2Num& num) {
  is.precision(std::numeric_limits<double>::max_digits10);
  std::string str;
  is >> str;
  if (str == "LOG2_ZERO") {
    num.is_log_zero_ = true;
  } else {
    num.is_log_zero_ = false;
    num.exponent_ = std::stod(str.substr(2));
  }
  return is;
}
