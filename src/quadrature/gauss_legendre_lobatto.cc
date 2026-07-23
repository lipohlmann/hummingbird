#include "quadrature/gauss_legendre_lobatto.h"

#include "math/legendre_polynomials.h"

namespace hummingbird::quadrature {
GaussLegendreLobatto::GaussLegendreLobatto(const size_t n_points)
    : QuadratureBase<double>(ComputeAbscissas(n_points)) {}

std::vector<double> GaussLegendreLobatto::ComputeAbscissas(
    const size_t n_points) {
  std::vector<double> abscissas(n_points);

  abscissas[0] = -1.0;
  abscissas.back() = 1.0;

  auto legendre_derivative_roots = math::AllLegendrePrimeRoots(n_points - 1);
  for (auto i = 1; i < abscissas.size() - 1; i++)
    abscissas[i] = legendre_derivative_roots[i - 1];
  return abscissas;
}
}  // namespace hummingbird::quadrature
