// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_H_
#define HUMMINGBIRD_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_H_

#include "quadrature/angular/gauss_legendre_chebyshev.h"
#include "quadrature/angular/gauss_legendre_trapezoid.h"
#include "quadrature/quadrature_base.h"
#include "utils/enums.h"

namespace hummingbird {
/**
 * @brief Wrapper class that owns an angular QuadratureBase<Ordinate> of the
 * type corresponding to the given AngularQuadSet
 *
 */
class AngularQuadratureSet {
 public:
  /**
   * @brief Construct a new AngularQuadratureSet object, creating the
   * concrete quadrature set matching quad_set
   *
   * @param quad_set Angular quadrature set to use
   * @param n_azim Number of azimuthal quadrature points
   * @param n_polar Number of polar quadrature points
   * @throw std::runtime_error if n_azim is not strictly positive, n_polar is
   * less than 2, n_polar is odd, or quad_set is not supported
   */
  AngularQuadratureSet(const AngularQuadSet quad_set, const size_t n_azim,
                       const size_t n_polar);

  /**
   * @brief Get the underlying quadrature set object
   *
   * @return const QuadratureBase<Ordinate>*
   */
  const QuadratureBase<Ordinate>* get() const {
    return angular_quadrature_set_.get();
  }

 private:
  /// @brief Concrete angular quadrature set for the chosen AngularQuadSet
  std::unique_ptr<QuadratureBase<Ordinate>> angular_quadrature_set_;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_H_
