#include "quadrature/gauss_lobatto_legendre.h"

#include "math/legendre_polynomials.h"

namespace hummingbird {
GaussLobattoLegendre::GaussLobattoLegendre(const size_t n_points)
    : QuadratureBase<double>(ComputeAbscissas(n_points)) {
  CreateWeightMap();
}

std::vector<double> GaussLobattoLegendre::ComputeAbscissas(
    const size_t n_points) {
  std::vector<double> abscissas(n_points);

  abscissas[0] = -1.0;
  abscissas.back() = 1.0;

  auto legendre_derivative_roots = AllLegendrePrimeRoots(n_points - 1);
  for (auto i = 1; i < abscissas.size() - 1; i++)
    abscissas[i] = legendre_derivative_roots[i - 1];
  return abscissas;
}

double GaussLobattoLegendre::ComputeWeight(const size_t k, const size_t n) {
  double double_n = static_cast<double>(n);
  double leading_coefficient = 2.0 / (double_n * (double_n - 1.0));
  double legendre_weight = LegendrePolynomial(n - 1, abscissas_.at(k));
  double weight = leading_coefficient / legendre_weight / legendre_weight;
  return weight;
}
}  // namespace hummingbird
