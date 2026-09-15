// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_QUADRATURE_GAUSS_LOBATTO_LEGENDRE_H_
#define HUMMINGBIRD_QUADRATURE_GAUSS_LOBATTO_LEGENDRE_H_

#include "math/legendre_polynomials.h"
#include "quadrature/quadrature_base.h"
namespace hummingbird {

/**
 * @brief Class defining a 1D Gauss-Legendre-Lobatto quadrature set on [-1,1].
 * The quadrature set approximates the integral using a set of weights
 * (\f$w_k\f$) and abscissas (\f$\xi_k\f$):
 *
 * \f[
 * \int_{-1}^{1}u(x)dx \approx \sum_{i=0}^{N-1}w_k u(\xi_k).
 * \f]
 *
 * Formulae for this class were taken from the textbook "High-Order Methods for
 * Incompressible Fluid Flow" by Deville, Fischer, and Mund.
 * https://doi.org/10.1017/CBO9780511546792
 *
 */
class GaussLobattoLegendre : public QuadratureBase<double> {
 public:
  GaussLobattoLegendre(const size_t n_points);

  /**
   * @brief Get the derivative of the specific Lagrange polynomial at the
   * specified node. Note that both the node index and the polynomial index are
   * zero-indexed.
   *
   * @param node_idx GLL node index
   * @param polynomial_idx Lagrange polynomial index
   * @return double
   */
  double GetLagrangeDerivative(const size_t node_idx,
                               const size_t polynomial_idx) const {
    size_t flattened_idx = node_idx * this->n_points() + polynomial_idx;
    return lagrange_derivatives_.at(flattened_idx);
  };

 private:
  /// @brief Derivatives of the Lagrange polynomials at the GLL nodes. These are
  /// stored as a flattened array and are indexed using the node number
  /// and the polynomial number in GetLagrangeDerivative. See the extended
  /// description in ComputeLagrangeDerivatives().
  std::vector<double> lagrange_derivatives_;

  /**
   * @brief Computes and sets the Lagrange polynomial derivatives at the GLL
   * nodes. The vector lagrange_derivatives_ is constructed with the
   * polynomial index being the "fast" counting index, and the node number
   * being the "slow" counting index. That is, the derivative of the N=1
   * polynomial at the 3rd GLL node of 5 has a flattened index of 16.
   *
   */
  void ComputeLagrangeDerivatives();

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
}  // namespace hummingbird

#endif  // HUMMINGBIRD_QUADRATURE_GAUSS_LOBATTO_LEGENDRE_H_
