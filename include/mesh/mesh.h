// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_MESH_MESH_H_
#define HUMMINGBIRD_MESH_MESH_H_

#include <string>
#include <vector>

#include "mesh/element.h"
#include "mesh/node.h"
#include "quadrature/gauss_lobatto_legendre.h"

namespace hummingbird {

/**
 * @brief Class defing a mesh
 *
 */
class Mesh {
 public:
  /**
   * @brief Construct a new Mesh object
   *
   */
  Mesh() = default;

  /**
   * @brief Construct a new Mesh object from a gmsh .msh file
   *
   * @param msh_file gmsh .msh file
   */
  Mesh(const std::string& msh_file);

  /**
   * @brief Add node to Mesh
   *
   * @param node Node
   */
  void AddNode(const Node& node);

  /**
   * @brief Add nodes in a vector to Mesh
   *
   * @param nodes Nodes
   */
  void AddNodes(const std::vector<Node>& nodes);

  /**
   * @brief Add element to Mesh
   *
   * @param element Element
   */
  void AddElement(std::unique_ptr<Element> element);

  /**
   * @brief Create interior nodes on elements using Gauss-Lobatto-Legendre
   * quadrature set and add to Mesh
   *
   * @param gll_quadrature Gauss-Lobatto-Legendre quadrature set
   */
  void CreateInteriorElementNodes(const GaussLobattoLegendre& gll_quadrature);

 private:
  /// @brief Nodes in the mesh
  std::vector<Node> nodes_;

  /// @brief Elements in the mesh
  std::vector<std::unique_ptr<Element>> elements_;

  void ReadGMSH(const std::string& msh_file);

  void BuildSegments();

  /**
   * @brief Renumber nodes in mesh to keep node IDs near each other in a single
   * element
   *
   */
  void RenumberNodes();

  /**
   * @brief Check that all Node IDs are unique and continuous from 0 to N-1 for
   * N total nodes
   *
   * @throw std::runtime_error Prints expected ID, found ID, and previous ID.
   *
   */
  void CheckNodeIDs();
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_MESH_MESH_H_
