#ifndef HUMMINGBIRD_TESTS_QUADRATURE_HELPERS_H_
#define HUMMINGBIRD_TESTS_QUADRATURE_HELPERS_H_

#include <functional>
#include <vector>

#include "quadrature/quadrature_base.h"

namespace hummingbird::quadrature {

template <typename Q>
std::vector<QuadraturePair> EvaluateAt(const Q& quad,
                                       const std::function<double(double)>& f) {
  std::vector<QuadraturePair> pairs;
  pairs.reserve(quad.abscissas().size());
  for (auto i = 0; i < quad.abscissas().size(); i++) {
    pairs.push_back(QuadraturePair{i, f(quad.GetAbscissa(i))});
  }
  return pairs;
}

double AnalyticMonomialIntegral(const unsigned int k);
}  // namespace hummingbird::quadrature

#endif  // HUMMINGBIRD_TESTS_QUADRATURE_HELPERS_H_