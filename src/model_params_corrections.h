#include "model_params_corrections.h"

#include "fast5/src/fast5.hpp"

using ::fast5::Model_Entry;
using ::fast5::Model_Parameters;

struct Gaussian {
  double mu_;
  double sigma_;
};

struct Gamma {
  double mean_;
  double stdv_;
};

Gaussian scaleGaussianCurrentLevel(const Gaussian& gaussian,
                                   const Model_Parameters& params) {
  double mu = gaussian.mu_ * params.scale + params.shift;
  double sigma = gaussian.sigma_ * params.var;

  return {mu, sigma};
}

Gamma scaleGammaNoiseLevel(const Gamma& gamma, const Model_Parameters& params) {
  double mean = gamma.mean_ * params.scale_sd;
  double stdv = gamma.stdv_ * params.var_sd;

  return {mean, stdv};
}
