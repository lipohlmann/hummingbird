#include "../tests/quadrature/helpers.h"

namespace hummingbird::quadrature {
std::vector<QuadraturePair> GLEvaluateAt(
    const GaussLegendre& quad, const std::function<double(double)>& f) {
  std::vector<QuadraturePair> pairs;
  pairs.reserve(quad.abscissas().size());
  for (auto i = 0; i < quad.abscissas().size(); i++) {
    pairs.push_back(QuadraturePair{i, f(quad.GetAbscissa(i))});
  }
  return pairs;
}

std::vector<QuadraturePair> GLLEvaluateAt(
    const GaussLegendreLobatto& quad, const std::function<double(double)>& f) {
  std::vector<QuadraturePair> pairs;
  pairs.reserve(quad.abscissas().size());
  for (auto i = 0; i < quad.abscissas().size(); i++) {
    pairs.push_back(QuadraturePair{i, f(quad.GetAbscissa(i))});
  }
  return pairs;
}
// Closed-form value of \int_{-1}^{1} x^k dx.
double AnalyticMonomialIntegral(const unsigned int k) {
  if (k % 2 == 1) {
    return 0.0;
  }
  return 2.0 / static_cast<double>(k + 1);
}
}  // namespace hummingbird::quadrature
