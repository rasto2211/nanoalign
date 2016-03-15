#include "log2_num.h"

Log2Num::Log2Num(double val) {
  if (val == 0) {
    is_log_zero_ = true;
  } else {
    exponent_ = log2(val);
    is_log_zero_ = false;
  }
}

double Log2Num::value() const {
  if (is_log_zero_) return 0;
  return exp2(exponent_);
}

void Log2Num::setExponent(double exponent) {
  is_log_zero_ = false;
  exponent_ = exponent;
}

Log2Num Log2Num::operator*(const Log2Num& num) const {
  Log2Num res = *this;
  res *= num;
  return res;
}

Log2Num& Log2Num::operator*=(const Log2Num& num) {
  if (num.isLogZero() || this->isLogZero()) {
    this->is_log_zero_ = true;
  } else {
    this->exponent_ += num.exponent_;
  }

  return *this;
}

Log2Num Log2Num::operator+(const Log2Num& num) const {
  Log2Num res = *this;
  res += num;
  return res;
}

Log2Num& Log2Num::operator+=(const Log2Num& num) {
  if (this->isLogZero()) {
    *this = num;
  } else if (!num.isLogZero()) {
    if (this->exponent_ > num.exponent_) {
      this->setExponent(this->exponent_ +
                        log2(1 + exp2(num.exponent_ - this->exponent_)));
    } else {
      this->setExponent(this->exponent_ +
                        log2(1 + exp2(this->exponent_ - num.exponent_)));
    }
  }

  return *this;
}

bool Log2Num::operator<(const Log2Num& num) const {
  if (is_log_zero_ && !num.is_log_zero_) return true;
  if (num.is_log_zero_) return false;
  return exponent_ < num.exponent_;
}

bool Log2Num::operator>(const Log2Num& num) const {
  if (num.is_log_zero_ && !is_log_zero_) return true;
  if (is_log_zero_) return false;
  return exponent_ > num.exponent_;
}
