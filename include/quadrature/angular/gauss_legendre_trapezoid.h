// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_TRAPEZOID_H_
#define HUMMINGBIRD_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_TRAPEZOID_H_

#include "quadrature/angular/ordinate.h"
#include "quadrature/gauss_legendre.h"
#include "quadrature/quadrature_base.h"
namespace hummingbird {

/**
 * @brief Class defines a Gauss-Legendre-Trapezoid rule on the unit sphere,
 * where a Gauss-Legendre rule is used on the polar coordinate, and a trapezoid
 * rule is used on the azimuthal.
 *
 */
class GaussLegendreTrapezoid : public QuadratureBase<Ordinate> {
 public:
  /**
   * @brief Construct a new Gauss Legendre Trapezoid object
   *
   * @param n_azim Number of azimuthal quadrature points
   * @param n_polar Number of polar quadrature points
   * @param polar_measured_from_x If true, each ordinate is built so its
   * polar Gauss-Legendre root ends up in Ordinate::x() (via
   * Ordinate(acos(root), pi/2)) instead of Ordinate::z(). Used for 1D
   * problems along x, where x() (not z()) is the streaming/BC direction
   * cosine every consumer (Segment, ParsedVolumetricSource, CGProblem)
   * reads; azim is otherwise degenerate for n_azim=1, since sin(polar) is
   * the same for a Gauss-Legendre root and its negative. Defaults to false
   * (the standard sphere convention, needed for future 2D/3D use).
   */
  GaussLegendreTrapezoid(const size_t n_azim, const size_t n_polar,
                        const bool polar_measured_from_x = false);
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_TRAPEZOID_H_
