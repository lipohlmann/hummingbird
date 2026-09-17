#ifndef HUMMINGBIRD_TESTS_QUADRATURE_HELPERS_H_
#define HUMMINGBIRD_TESTS_QUADRATURE_HELPERS_H_

#include <functional>
#include <vector>

#include "quadrature/quadrature_base.h"

namespace hummingbird {

/**
 * @brief Evaluate a function at each abscissa of a quadrature set, returning
 * the QuadraturePair objects needed to integrate it
 *
 * @tparam Q Quadrature set type
 * @param quad Quadrature set
 * @param f Function to evaluate at each abscissa
 * @return std::vector<QuadraturePair>
 */
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

/**
 * @brief Compute the analytic value of \f$\int_{-1}^1 x^k dx\f$, used to
 * check quadrature sets against monomial test functions
 *
 * @param k Monomial order
 * @return double
 */
double AnalyticMonomialIntegral(const unsigned int k);
}  // namespace hummingbird

#endif  // HUMMINGBIRD_TESTS_QUADRATURE_HELPERS_H_