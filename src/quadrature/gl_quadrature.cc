#include "quadrature/gl_quadrature.h"

#include "math/legendre_polynomials.h"

namespace hummingbird::quadrature {
GLQuadrature::GLQuadrature(const unsigned int n)
    : Quadrature<double>(math::AllLegendreRoots(n)) {}
}  // namespace hummingbird::quadrature
