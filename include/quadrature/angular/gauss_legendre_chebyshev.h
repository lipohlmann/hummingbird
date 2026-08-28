// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_CHEBYSHEV_H_
#define HUMMINGBIRD_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_CHEBYSHEV_H_

#include "quadrature/angular/ordinate.h"
#include "quadrature/gauss_chebyshev.h"
#include "quadrature/gauss_legendre.h"
#include "quadrature/quadrature_base.h"

namespace hummingbird::quadrature::angular {
/**
 * @brief Class defines a Gauss-Legendre-Chebyshev quadrature set on the unit
 * sphere, where a Gauss-Legendre quadrature set is used along the polar
 * direction, and a Gauss-Chebyshev quadrature set is used along the azimuthal
 * direction.
 *
 */
class GaussLegendreChebyshev : public quadrature::QuadratureBase<Ordinate> {
 public:
  /**
   * @brief Construct a new Gauss Legendre Chebyshev object
   *
   * @param n_half_azim Number of directions in a half-circle (directions will
   * be mirrored, and thus doubled)
   * @param n_polar Number of directions on the polar axis
   */
  GaussLegendreChebyshev(const size_t n_half_azim, const size_t n_polar);

 private:
};

}  // namespace hummingbird::quadrature::angular

#endif  // HUMMINGBIRD_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_CHEBYSHEV_H_
