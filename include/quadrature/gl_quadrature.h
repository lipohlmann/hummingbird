#ifndef HUMMINGBIRD_QUADRATURE_GL_QUADRATURE_H
#define HUMMINGBIRD_QUADRATURE_GL_QUADRATURE_H

#include "quadrature/quadrature.h"
namespace hummingbird::quadrature {

/**
 * @brief Defines a Gauss-Legendre quadrature set on the interval [-1,1]
 *
 */
class GaussLegendreQuadrature : public Quadrature<double> {
 public:
  /**
   * @brief Construct a new GaussLegendreQuadrature object
   *
   * @param n_points Number of quadrature points to be created
   */
  GaussLegendreQuadrature(const unsigned int n_points);

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
