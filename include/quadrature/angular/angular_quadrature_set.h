// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_H_
#define HUMMINGBIRD_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_H_

#include "quadrature/angular/gauss_legendre_chebyshev.h"
#include "quadrature/angular/gauss_legendre_trapezoid.h"
#include "quadrature/quadrature_base.h"
#include "utils/enums.h"

namespace hummingbird {
class AngularQuadratureSet {
 public:
  AngularQuadratureSet(const AngularQuadSet quad_set, const size_t n_azim,
                       const size_t n_polar);

  const QuadratureBase<Ordinate>* get() const {
    return angular_quadrature_set_.get();
  }

 private:
  std::unique_ptr<QuadratureBase<Ordinate>> angular_quadrature_set_;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_SET_H_
