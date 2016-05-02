#include <cassert>

#include "log2_num.h"

#include <glog/logging.h>

const double kEpsilon = 1.0e-15;

Log2Num::Log2Num(double val) { exponent_ = log2(val); }

double Log2Num::value() const {
  if (isLogZero()) return 0;
  return exp2(exponent_);
}

void Log2Num::setExponent(double exponent) { exponent_ = exponent; }

Log2Num Log2Num::operator*(const Log2Num& num) const {
  Log2Num res = *this;
  res *= num;
  return res;
}

Log2Num& Log2Num::operator*=(const Log2Num& num) {
  if (num.isLogZero() || this->isLogZero()) {
    this->exponent_ = -HUGE_VAL;
  } else {
    this->exponent_ += num.exponent_;
  }

  return *this;
}

Log2Num Log2Num::operator+(const Log2Num& num) const {
  Log2Num res(*this);
  res += num;
  return res;
}

Log2Num& Log2Num::operator+=(const Log2Num& num) {
  if (this->isLogZero()) {
    *this = num;
  } else if (!num.isLogZero()) {
    if (this->exponent_ > num.exponent_) {
      this->exponent_ += log2(1 + exp2(num.exponent_ - this->exponent_));
    } else {
      this->exponent_ =
          num.exponent_ + log2(1 + exp2(this->exponent_ - num.exponent_));
    }
  }

  return *this;
}

bool Log2Num::operator<(const Log2Num& num) const {
  if (isLogZero() && !num.isLogZero()) return true;
  if (num.isLogZero()) return false;
  return exponent_ < num.exponent_;
}

bool Log2Num::operator>(const Log2Num& num) const {
  if (num.isLogZero() && !isLogZero()) return true;
  if (isLogZero()) return false;
  return exponent_ > num.exponent_;
}

bool Log2Num::operator==(const Log2Num& num) const {
  double this_val = this->value();
  double num_val = num.value();
  return fabs(this_val - num_val) < kEpsilon;
}

bool Log2Num::operator!=(const Log2Num& num) const { return !(*this == num); }

Log2Num& Log2Num::operator/=(const Log2Num& num) {
  assert(!num.isLogZero());

  if (!this->isLogZero()) this->exponent_ -= num.exponent_;

  return *this;
}

Log2Num Log2Num::operator/(const Log2Num& num) const {
  Log2Num res = *this;
  res /= num;
  return res;
}
