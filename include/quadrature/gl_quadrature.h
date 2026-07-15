#ifndef HUMMINGBIRD_QUADRATURE_GL_QUADRATURE_H
#define HUMMINGBIRD_QUADRATURE_GL_QUADRATURE_H

#include "quadrature/quadrature.h"
namespace hummingbird::quadrature {
class GLQuadrature : public Quadrature<double> {
 public:
  GLQuadrature(const unsigned int n_ordinates);

 private:
  /**
   * @brief Compute the weight for abscissa value k for a n-th order
   * Gauss-Legendre polynomial
   *
   * @param k Abscissa number
   * @param n Order of Gauss-Legendre polynomial
   * @return Weight corresponding to the k-th abscissa values
   */
  double ComputeWeight(const auto k, const auto n);

  /**
   * @brief Create the weight map
   *
   */
  void CreateWeightMap();
};
}  // namespace hummingbird::quadrature

#endif  // HUMMINGBIRD_QUADRATURE_GL_QUADRATURE_H
