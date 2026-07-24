#include "quadrature/gauss_chebyshev.h"

#include <cmath>
namespace hummingbird::quadrature {
GaussChebyshev::GaussChebyshev(const size_t n_points)
    : QuadratureBase<double>(ComputeChebyshevAbscissas(n_points)) {
  CreateWeightMap();
}

std::vector<double> GaussChebyshev::ComputeChebyshevAbscissas(
    const size_t n_points) {
  std::vector<double> abscissas(n_points);
  for (auto i = 0; i < n_points; i++)
    abscissas[i] = std::cos((2.0 * static_cast<double>(i) + 1.0) /
                            (2.0 * static_cast<double>(n_points)) * M_PI);
  return abscissas;
}

double GaussChebyshev::ComputeWeight(const size_t k, const size_t n) {
  return M_PI / static_cast<double>(n);
}
}  // namespace hummingbird::quadrature
