#include "quadrature/angular/angular_quadrature_set.h"

#include <stdexcept>

namespace hummingbird {
AngularQuadratureSet::AngularQuadratureSet(const AngularQuadSet quad_set,
                                           const size_t n_azim,
                                           const size_t n_polar) {
  if (n_azim < 1)
    throw std::runtime_error(
        "Number of azimuthal angles must be strictly positive");
  if (n_polar < 2)
    throw std::runtime_error("Number of polar angles must be at least 2.");
  if (n_polar % 2 != 0)
    throw std::runtime_error("Number of polar angles must be even.");

  switch (quad_set) {
    case AngularQuadSet::GL:
      angular_quadrature_set_ =
          std::make_unique<GaussLegendreTrapezoid>(1, n_polar);
      break;
    case AngularQuadSet::GLT:
      angular_quadrature_set_ =
          std::make_unique<GaussLegendreTrapezoid>(n_azim, n_polar);
      break;
    case AngularQuadSet::GLC:
      angular_quadrature_set_ =
          std::make_unique<GaussLegendreChebyshev>(n_azim, n_polar);
      break;
    default:
      throw std::runtime_error(
          "AngularQuadSet passed to AngularQuadratureSet constructor not "
          "supported.");
      break;
  }
}
}  // namespace hummingbird
