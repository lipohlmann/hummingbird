#include "quadrature/gauss_lobatto_legendre.h"

#include <cassert>
#include <stdexcept>

#include "math/legendre_polynomials.h"

namespace hummingbird {
GaussLobattoLegendre::GaussLobattoLegendre(const size_t n_points)
    : QuadratureBase<double>(ComputeAbscissas(n_points)) {
  CreateWeightMap();
  ComputeLagrangeDerivatives();
}

std::vector<double> GaussLobattoLegendre::ComputeAbscissas(
    const size_t n_points) {
  if (n_points < 2)
    throw std::runtime_error(
        "n_points in GaussLobattoLegendre must be at least 2.");
  std::vector<double> abscissas(n_points);

  abscissas[0] = -1.0;
  abscissas.back() = 1.0;

  auto legendre_derivative_roots = AllLegendrePrimeRoots(n_points - 1);
  for (auto i = 1; i < abscissas.size() - 1; i++)
    abscissas[i] = legendre_derivative_roots[i - 1];
  return abscissas;
}

void GaussLobattoLegendre::ComputeLagrangeDerivatives() {
  auto n_points = this->n_points();
  lagrange_derivatives_.reserve(n_points * n_points);
  for (auto i = 0; i < n_points; i++) {
    for (auto k = 0; k < n_points; k++) {
      if (i == 0 && k == 0)
        lagrange_derivatives_.push_back(-(n_points * (n_points - 1) / 4.0));
      else if (i == (n_points - 1) && k == (n_points - 1))
        lagrange_derivatives_.push_back(n_points * (n_points - 1) / 4.0);
      else if (i != k) {
        auto abscissa_i = abscissas_.at(i);
        auto abscissa_k = abscissas_.at(k);
        double legendre_ratio = LegendrePolynomial(n_points - 1, abscissa_i) /
                                LegendrePolynomial(n_points - 1, abscissa_k);
        double abscissa_difference = abscissa_i - abscissa_k;
        lagrange_derivatives_.push_back(legendre_ratio / abscissa_difference);
      } else {
        lagrange_derivatives_.push_back(0.0);
      }
    }
  }
}

double GaussLobattoLegendre::ComputeWeight(const size_t k, const size_t n) {
  double double_n = static_cast<double>(n);
  double leading_coefficient = 2.0 / (double_n * (double_n - 1.0));
  double legendre_weight = LegendrePolynomial(n - 1, abscissas_.at(k));
  double weight = leading_coefficient / legendre_weight / legendre_weight;
  return weight;
}

double GaussLobattoLegendre::IntegrateGridFunction(
    const std::vector<double>& grid_function_vals) {
  assert(grid_function_vals.size() == abscissas_.size());
  double sum = 0.0;
  for (auto i = 0; i < grid_function_vals.size(); i++)
    sum += grid_function_vals[i] * weight_map_.at(i);
  return sum;
}
}  // namespace hummingbird
