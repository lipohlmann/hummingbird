// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_MESH_ELEMENT_H_
#define HUMMINGBIRD_MESH_ELEMENT_H_

#include <armadillo>
#include <memory>
#include <vector>

#include "banks/material_bank.h"
#include "node.h"
#include "physics/material.h"
#include "quadrature/gauss_lobatto_legendre.h"

namespace hummingbird {
/**
 * @brief Defines a subset of the domain (an "element")
 *
 */
class Element {
 public:
  /**
   * @brief Construct a new Element object
   *
   * @param material_id Material ID
   * @param source_id Source ID
   */
  Element(const int material_id, const int source_id);

  /**
   * @brief Create nodes by mapping the Gauss-Lobatto-Legendre quadrature set
   * from the reference to the real domain
   *
   * @param existing_nodes Vector of nodes that already exist in the mesh.
   * Starting ID for new nodes will be equal to the number of nodes that already
   * exist to ensure no duplicates
   * @param gll_quadrature GaussLobattoLegendre quadrature set
   */
  virtual std::vector<Node> CreateInteriorNodes(

      const std::vector<Node>& existing_nodes,
      const GaussLobattoLegendre& gll_quadrature) = 0;

  /**
   * @brief Get node IDs
   *
   * @return std::vector<size_t>
   */
  std::vector<size_t> node_ids() const { return node_ids_; }

  /**
   * @brief Get Material ID
   *
   * @return int
   */
  int material_id() const { return material_id_; }

  /**
   * @brief Get Source ID
   *
   * @return int
   */
  int source_id() const { return source_id_; }

  /**
   * @brief Set the new Node ID
   *
   * @param prev_id Previous ID (currently stored in object)
   * @param new_id New ID
   */
  void SetNewNodeID(const size_t prev_id, const size_t new_id);

  //   virtual arma::Mat<double> LocalStiffnessMatrix() = 0;

  /**
   * @brief Construct the sparse *diagonal* local mass matrix using spectral
   * elements on a Gauss-Lobatto-Legendre grid, defined as:
   *
   * \f[
   * \Sigma_t^e
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
   * @param material_bank Material bank
   * @return arma::SpMat<double>
   */
  virtual arma::SpMat<double> LocalMassMatrix(
      const GaussLobattoLegendre& gll_quad,
      const MaterialBank& material_bank) = 0;

 protected:
  /// @brief Material ID
  const int material_id_;

  /// @brief Source ID
  const int source_id_;

  /// @brief IDs of nodes defining the element
  std::vector<size_t> node_ids_;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_MESH_ELEMENT_H_
