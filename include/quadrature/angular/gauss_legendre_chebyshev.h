#ifndef HUMMINGBIRD_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_CHEBYSHEV_H_
#define HUMMINGBIRD_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_CHEBYSHEV_H_

#include "quadrature/angular/ordinate.h"
#include "quadrature/gauss_chebyshev.h"
#include "quadrature/gauss_legendre.h"
#include "quadrature/quadrature_base.h"

namespace hummingbird::quadrature::angular {
class GaussLegendreChebyshev : public quadrature::QuadratureBase<Ordinate> {
 public:
  GaussLegendreChebyshev(const size_t n_half_azim, const size_t n_polar);

 private:
};

}  // namespace hummingbird::quadrature::angular

#endif  // HUMMINGBIRD_QUADRATURE_ANGULAR_GAUSS_LEGENDRE_CHEBYSHEV_H_
