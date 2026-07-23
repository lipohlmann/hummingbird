#ifndef HUMMINGBIRD_QUADRATURE_GAUSS_LEGENDRE_LOBATTO_H_
#define HUMMINGBIRD_QUADRATURE_GAUSS_LEGENDRE_LOBATTO_H_

#include "quadrature/quadrature_base.h"

namespace hummingbird::quadrature {

/**
 * @brief Class defining a 1D Gauss-Legendre-Lobatto quadrature set on [-1,1].
 * The quadrature set approximates the integral using a set of weights
 * (\f$w_k\f$) and abscissas (\f$\xi_k\f$):
 *
 * \f[
 * \int_{-1}^{1}u(x)dx \approx \sum_{i=0}^{N-1}w_k u(\xi_k).
 * \f]
 *
 * Formulae for this class was taken from the textbook "High-Order Methods for
 * Incompressible Fluid Flow" by Deville, Fischer, and Mund.
 * https://doi.org/10.1017/CBO9780511546792
 *
 */
class GaussLegendreLobatto : public QuadratureBase<double> {
 public:
  GaussLegendreLobatto(const size_t n_points);

 private:
  /**
   * @brief Computes the abscissa values, given by:
   *
   * \f[
   * \xi_k = \begin{cases}
   * -1, & k=0\\
   * \text{zeros of } P'_{N-1}, & 1\leq k\leq N-2\\
   * 1,& k=N-1
   * \end{cases}
   * \f]
   *
   * @param n_points Total number of abscissa points
   * @return std::vector<double>
   */
  std::vector<double> ComputeAbscissas(const size_t n_points);

  /**
   * @brief Computes the weights for a Gauss-Legendre-Lobatto quadrature scheme
   * using:
   *
   * \f[
   * w_k=\frac{2}{(N-1)N}\frac{1}{[P_{N-1}(\xi_k)]^2}
   * \f]
   *
   * @param k Abscissa index
   * @param n Total number of points
   * @return Quadrature weight associated with abscissa k
   */
  double ComputeWeight(const size_t k, const size_t n) override;
};
}  // namespace hummingbird::quadrature

#endif  // HUMMINGBIRD_QUADRATURE_GAUSS_LEGENDRE_LOBATTO_H_
