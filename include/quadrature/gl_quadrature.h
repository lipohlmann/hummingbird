#ifndef HUMMINGBIRD_QUADRATURE_GL_QUADRATURE_H
#define HUMMINGBIRD_QUADRATURE_GL_QUADRATURE_H

#include "quadrature/quadrature.h"
namespace hummingbird::quadrature {
class GLQuadrature : public Quadrature<double> {
 public:
  GLQuadrature(const unsigned int n_ordinates);
};
}  // namespace hummingbird::quadrature

#endif  // HUMMINGBIRD_QUADRATURE_GL_ANGULAR_QUADRATURE_H
