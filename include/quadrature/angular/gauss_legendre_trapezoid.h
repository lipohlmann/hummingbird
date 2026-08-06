#ifndef HUMMINGBIRD_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_TRAPEZOID_H_
#define HUMMINGBIRD_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_TRAPEZOID_H_

#include "quadrature/angular/ordinate.h"
#include "quadrature/gauss_legendre.h"
#include "quadrature/quadrature_base.h"
namespace hummingbird::quadrature::angular {

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
   */
  GaussLegendreTrapezoid(const size_t n_azim, const size_t n_polar);
};
}  // namespace hummingbird::quadrature::angular

#endif  // HUMMINGBIRD_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_TRAPEZOID_H_
