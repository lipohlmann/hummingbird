#include "quadrature/angular/gauss_legendre_trapezoid.h"

#include <cmath>
namespace hummingbird::quadrature::angular {
GaussLegendreTrapezoid::GaussLegendreTrapezoid(const size_t n_azim,
                                               const size_t n_polar) {
  abscissas_.reserve(n_azim * n_polar);
  GaussLegendre gl_polar_quad(n_polar);

  double delta_azim = 2.0 * M_PI / (static_cast<double>(n_azim - 1));
  size_t counter = 0;
  for (auto i = 0; i < n_polar; i++) {
    double polar = gl_polar_quad.GetAbscissa(i);
    for (double j = 0; j < n_azim; j++) {
      double azim = j * delta_azim;

      Ordinate new_ordinate(std::acos(azim), std::acos(polar));
      abscissas_.push_back(std::move(new_ordinate));

      double weight = gl_polar_quad.GetWeight(i) * delta_azim;
      weight_map_.insert({counter, weight});
      counter++;
    }
  }
}
}  // namespace hummingbird::quadrature::angular
