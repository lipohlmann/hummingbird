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
#include "quadrature/angular/ordinate.h"
#include "quadrature/gauss_lobatto_legendre.h"

namespace hummingbird {
class Mesh;

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

  /**
   * @brief Replace all node IDs at once, e.g. when renumbering. Unlike
   * SetNewNodeID, this does not search by value, so it is safe to use when
   * remapped IDs could otherwise collide with not-yet-remapped IDs.
   *
   * @param node_ids New node IDs, in the same order as the existing ones
   */
  void SetNodeIDs(std::vector<size_t> node_ids) {
    node_ids_ = std::move(node_ids);
  }

  virtual arma::Col<double> LocalForcingVector(
      const GaussLobattoLegendre& gll_quad, const Mesh& mesh,
      const size_t ordinate_index) = 0;

  virtual arma::Mat<double> LocalStiffnessMatrix(
      const GaussLobattoLegendre& gll_quad, const MaterialBank& material_bank,
      const Ordinate& ordinate) = 0;

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
