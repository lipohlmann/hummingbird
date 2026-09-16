// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_MESH_SEGMENT_H_
#define HUMMINGBIRD_MESH_SEGMENT_H_

#include <array>

#include "banks/material_bank.h"
#include "mesh/element.h"
#include "mesh/mesh.h"
#include "mesh/node.h"

namespace hummingbird {
class Segment : public Element {
 public:
  /**
   * @brief Construct a new Segment object
   *
   * @param boundary_node_ids Node IDs defining the Segment bounds
   * @param material_id Material ID
   * @param source_id Source ID
   */
  Segment(const std::array<size_t, 2> boundary_node_ids, const int material_id,
          const int source_id, const Mesh& mesh);

  /**
   * @brief Create interior Node objects using the Gauss-Lobatto-Legendre
   * quadrature set. The mapping from the reference to the real domain was based
   * on the equations provided by the book Computational Seismology by Heiner
   * Igel
   * (https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://www.geokniga.org/bookfiles/geokniga-computationalseismologyapracticalintroductionbyigelheinerz-liborg.pdf&ved=2ahUKEwiMuanq7_OVAxXVN4YAHVgAD1QQFnoECBUQAQ&usg=AOvVaw35g1Tl3ZSFUMusz8LKQxHk).
   * For a line segment defined by two endpoints, the mapping is:
   *
   * \f[
   * \mathbf{r}(\xi) = \mathbf{r}_1+\left( \frac{\xi+1}{2} \right)\left(
   * \mathbf{r}_2-\mathbf{r}_1 \right) \f]
   *
   * where \f$x\in[-1,1]\f$.
   *
   * @param existing_nodes Nodes already existing in Mesh
   * @param gll_quadrature Gauss-Lobatto-Legendre quadrature set
   * @return Additional nodes to be defined along the GLL quadrature set. Will
   * be of length N-2 (the two endpoint account for the remaining two nodes)
   */
  std::vector<Node> CreateInteriorNodes(
      const std::vector<Node>& existing_nodes,
      const GaussLobattoLegendre& gll_quadrature) override;

  /**
   * @brief Construct the local forcing vector, which is:
   *
   * \f[
   * f_i = \frac{h_e}{2}\int_{-1}^1Q_n(\xi)\ell_i (\xi)d\xi +
   * int_{-1}^1Q_n(\xi)\frac{d\ell_i}{d\xi}d\xi
   * \f]
   * which the integrals are again approximated using the GLL quadrature set:
   * \f[
   * f_i \approx \frac{h_e}{2} \sum_{k=0}^{N_x}\rho_k Q_n(\xi_k)\left[ \ell_i
   * (\xi_k) +\frac{d\ell_i}{d\xi}(\xi_k)\right]
   * \f]
   *
   * @param gll_quad
   * @param mesh
   * @param ordinate_index
   * @return arma::Col<double>
   */
  arma::Col<double> LocalForcingVector(const GaussLobattoLegendre& gll_quad,
                                       const Mesh& mesh,
                                       const size_t ordinate_index) override;

  /**
   * @brief Construct the dense local stiffness matrix using spectral elements
   * on a Gauss-Lobatto-Legendre grid, defined as:
   *
   * \f[
   * K_{ij}=\frac{\hat{\Omega}_x^2}{\Sigma_t^e}\frac{2}{h_e}\sum_{k=0}^{N_x}\rho_k
   * \frac{d\ell_i}{d\xi}(\xi_k)\frac{d\ell_j}{d\xi}(\xi_k)
   * \f]
   *
   * @return arma::Mat<double>
   * @todo See about cutting the loops in half by leveraging symmetry.
   */
  arma::Mat<double> LocalStiffnessMatrix(const GaussLobattoLegendre& gll_quad,
                                         const MaterialBank& material_bank,
                                         const Ordinate& ordinate) override;

  /**
   * @brief Construct the sparse *diagonal* local mass matrix using spectral
   * elements on a Gauss-Lobatto-Legendre grid, defined as:
   *
   * \f[
   * M_{ij}=\Sigma_t^e\frac{h_e}{2}\sum_{k=0}^{N_x}\rho_k\ell_i(\xi_k)\ell_j(\xi_k)
   * \f]
   * Which, using the cardinality of Lagrange polynomials, greatly simplifies
   * to:
   * \f[
   * \Sigma_t^eM_{ij}=\begin{cases}
   * 0, & i\neq j \\
   * \Sigma_t^e\frac{h_e}{2}\rho_i, &i=j
   * \end{cases}
   * \f]
   *
   * @param gll_quad GaussLobattoLegendre set
   * @param material_bank MaterialBank
   * @return arma::SpMat<double>
   */
  arma::SpMat<double> LocalMassMatrix(
      const GaussLobattoLegendre& gll_quad,
      const MaterialBank& material_bank) override;

  /**
   * @brief Get the spatial dimension of a Segment (always 1)
   *
   * @return unsigned int
   */
  unsigned int dimension() const override;

 private:
  std::array<size_t, 2> boundary_node_ids_;

  /// @brief Length of the element
  const double length_;

  /**
   * @brief Compute the length of the element on construction
   *
   * @param boundary_node_ids Boundary node IDs
   * @param mesh Mesh
   * @return double
   */
  double ComputeLength(const std::array<size_t, 2> boundary_node_ids,
                       const Mesh& mesh);
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_MESH_SEGMENT_H_
