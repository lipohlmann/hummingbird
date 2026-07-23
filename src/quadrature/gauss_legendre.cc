#include "quadrature/gauss_legendre.h"

#include "math/legendre_polynomials.h"

namespace hummingbird::quadrature {
GaussLegendre::GaussLegendre(const unsigned int n_points)
    : QuadratureBase<double>(math::AllLegendreRoots(n_points)) {
  CreateWeightMap();
}

double GaussLegendre::ComputeWeight(const size_t k, const size_t n) {
  double abscissa = abscissas_.at(k);
  double derivative_weight = math::LegendrePolynomialPrime(n, abscissa);
  return 2.0 / (1 - abscissa * abscissa) / derivative_weight /
         derivative_weight;
}
}  // namespace hummingbird::quadrature
