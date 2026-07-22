#include "quadrature/gl_quadrature.h"

#include "math/legendre_polynomials.h"

namespace hummingbird::quadrature {
GLQuadrature::GLQuadrature(const unsigned int n_points)
    : Quadrature<double>(math::AllLegendreRoots(n_points)) {
  CreateWeightMap();
}

void GLQuadrature::CreateWeightMap() {
  int order = abscissas_.size();
  for (auto i = 0; i < order; i++) {
    double weight = ComputeWeight(i, order);
    weight_map_.insert({i, weight});
  }
}

double GLQuadrature::ComputeWeight(const auto k, const auto n) {
  double abscissa = abscissas_.at(k);
  double derivative_weight = math::LegendrePolynomialPrime(n, abscissa);
  return 2.0 / (1 - abscissa * abscissa) / derivative_weight /
         derivative_weight;
}
}  // namespace hummingbird::quadrature
