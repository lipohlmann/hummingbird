#include "quadrature/angular/gauss_legendre_chebyshev.h"

#include <cmath>

#include "quadrature/gauss_chebyshev.h"
#include "quadrature/gauss_legendre.h"

namespace hummingbird::quadrature::angular {

GaussLegendreChebyshev::GaussLegendreChebyshev(const size_t n_half_azim,
                                               const size_t n_polar) {
  GaussChebyshev gc_azim_quad(n_half_azim);
  GaussLegendre gl_polar_quad(n_polar);
  size_t n_azim = 2 * n_half_azim;
  size_t n_total = n_polar * n_azim;
  abscissas_.reserve(n_total);

  size_t counter = 0;
  for (auto i = 0; i < n_polar; i++) {
    for (auto j = 0; j < n_half_azim; j++) {
      double azim = gc_azim_quad.GetAbscissa(j);
      double polar = gl_polar_quad.GetAbscissa(i);

      Ordinate positive(std::acos(azim), std::acos(polar));
      abscissas_.push_back(std::move(positive));
      double weight = gl_polar_quad.GetWeight(i) * gc_azim_quad.GetWeight(j);
      weight_map_.insert({counter, weight});
      counter++;

      Ordinate negative(std::acos(azim), std::acos(polar));
      abscissas_.push_back(std::move(negative));
      weight_map_.insert({counter, weight});
    }
  }
}

}  // namespace hummingbird::quadrature::angular
