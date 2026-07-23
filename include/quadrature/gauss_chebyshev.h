#ifndef HUMMINGBIRD_QUADRATURE_GAUSS_CHEBYSHEV_H_
#define HUMMINGBIRD_QUADRATURE_GAUSS_CHEBYSHEV_H_

#include "quadrature/quadrature_base.h"

namespace hummingbird::quadrature {
/**
 * @brief Class defines a Gauss-Chebyshev quadrature set, which approximates the
 * integral of the form:
 *
 * \f[
 * \int^{1}_{-1}\frac{u(x)}{\sqrt{1-x^2}}dx\approx \sum_{k=0}^{N-1}w_k u(\xi_k)
 * \f]
 *
 */
class GaussChebyshev : public QuadratureBase<double> {
 public:
  /**
   * @brief Construct a new Gauss Chebyshev object
   *
   * @param n_points Total number of points desired
   */
  GaussChebyshev(const size_t n_points);

 private:
  /**
   * @brief Computes the abscissa values for Gauss Chebyshev quadrature. These
   * are given by the explicit formula:
   *
   * \f[
   * \xi_k=\cos \left( \frac{(2k+1)\pi}{2N} \right)
   * \f]
   * for N total points.
   *
   * @param n_points
   * @return std::vector<double>
   */
  std::vector<double> ComputeChebyshevAbscissas(const size_t n_points);

  /**
   * @brief Computes the quadrature weight for abscissa k of n total points. In
   * Gauss-Chebyshev quadrature, all weights are equal and are given by:
   *
   * \f[
   * w_k=\frac{\pi}{n}
   * \f]
   *
   * @param k Abscissa number (not used in computation)
   * @param n Total number of points
   * @return Weight for abscissa k
   */
  double ComputeWeight(const size_t k, const size_t n) override;
};
}  // namespace hummingbird::quadrature

#endif  // HUMMINGBIRD_QUADRATURE_GAUSS_CHEBYSHEV_H_
