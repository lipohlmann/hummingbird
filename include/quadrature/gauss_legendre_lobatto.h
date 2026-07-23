#ifndef HUMMINGBIRD_QUADRATURE_GAUSS_LEGENDRE_LOBATTO_H_
#define HUMMINGBIRD_QUADRATURE_GAUSS_LEGENDRE_LOBATTO_H_

#include "quadrature/quadrature_base.h"

namespace hummingbird::quadrature {
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
  double ComputeWeight(const size_t k, const size_t n) override;
};
}  // namespace hummingbird::quadrature

#endif  // HUMMINGBIRD_QUADRATURE_GAUSS_LEGENDRE_LOBATTO_H_
