#ifndef HUMMINGBIRD_TESTS_QUADRATURE_HELPERS_H_
#define HUMMINGBIRD_TESTS_QUADRATURE_HELPERS_H_

#include <functional>
#include <vector>

#include "quadrature/gauss_legendre.h"
#include "quadrature/gauss_legendre_lobatto.h"
#include "quadrature/quadrature_base.h"

namespace hummingbird::quadrature {

std::vector<QuadraturePair> GLEvaluateAt(
    const GaussLegendre& quad, const std::function<double(double)>& f);

std::vector<QuadraturePair> GLLEvaluateAt(
    const GaussLegendreLobatto& quad, const std::function<double(double)>& f);

double AnalyticMonomialIntegral(const unsigned int k);
}  // namespace hummingbird::quadrature

#endif  // HUMMINGBIRD_TESTS_QUADRATURE_HELPERS_H_