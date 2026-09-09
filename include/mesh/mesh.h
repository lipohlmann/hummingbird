// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_MESH_MESH_H_
#define HUMMINGBIRD_MESH_MESH_H_

#include <fstream>
#include <string>
#include <unordered_map>
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

  /**
   * @brief Get the Nodes in the mesh
   *
   * @return const std::vector<Node>&
   */
  const std::vector<Node>& nodes() const { return nodes_; }

  /**
   * @brief Get the number of Elements in the mesh
   *
   * @return size_t
   */
  size_t n_elements() const { return elements_.size(); }

  /**
   * @brief Get an Element by index
   *
   * @param index Index of the element
   * @return const Element&
   */
  const Element& GetElement(const size_t index) const {
    return *elements_.at(index);
  }

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

  /**
   * @brief Intermediate state accumulated while reading a gmsh file, shared
   * across the section subroutines below.
   *
   */
  struct GmshReadState {
    /// @brief Map from Physical Group tag to its name (e.g. "material:foo")
    std::unordered_map<int, std::string> physical_names;

    /// @brief Map from point entity tag to its Physical Group tags
    std::unordered_map<int, std::vector<int>> point_physical_tags;

    /// @brief Map from curve entity tag to its Physical Group tags
    std::unordered_map<int, std::vector<int>> curve_physical_tags;

    /// @brief Map from gmsh node tag to the corresponding Node's ID in
    /// nodes_
    std::unordered_map<size_t, size_t> node_tag_to_id;
  };

  /**
   * @brief Construct a Mesh from a gmsh ASCII (format 4.1) .msh file. Reading
   * routes to a section-specific subroutine using an unordered_map keyed by
   * the gmsh section header (e.g. "$Nodes"). Only the sections needed to
   * build 1D meshes (Segment elements, from 2-node line elements and Physical
   * Point/Curve groups) are currently handled; support for 2D quad elements
   * can be added later by extending ReadElements without needing to
   * restructure this dispatch.
   *
   * @param msh_file Path to the gmsh .msh file
   */
  void ReadGMSH(const std::string& msh_file);

  /**
   * @brief Read the $PhysicalNames section into state.physical_names
   *
   * @param file gmsh file stream, positioned just after the section header
   * @param state Shared gmsh read state
   */
  void ReadPhysicalNames(std::ifstream& file, GmshReadState& state);

  /**
   * @brief Read the $Entities section into state.point_physical_tags and
   * state.curve_physical_tags
   *
   * @throw std::runtime_error if the mesh contains surface or volume
   * entities, since only 1D (point/curve) meshes are currently supported
   *
   * @param file gmsh file stream, positioned just after the section header
   * @param state Shared gmsh read state
   */
  void ReadEntities(std::ifstream& file, GmshReadState& state);

  /**
   * @brief Read the $Nodes section, adding Node objects to the mesh and
   * recording the gmsh tag to Node ID mapping in state.node_tag_to_id
   *
   * @param file gmsh file stream, positioned just after the section header
   * @param state Shared gmsh read state
   */
  void ReadNodes(std::ifstream& file, GmshReadState& state);

  /**
   * @brief Read the $Elements section, adding Segment elements to the mesh
   * for 2-node line elements. Point elements are consumed but otherwise
   * ignored, since they only mark boundary entities in gmsh.
   *
   * @throw std::runtime_error for any element type other than a point or a
   * 2-node line, since only 1D meshes are currently supported
   *
   * @param file gmsh file stream, positioned just after the section header
   * @param state Shared gmsh read state
   */
  void ReadElements(std::ifstream& file, GmshReadState& state);

  /**
   * @brief Find the material ID for a curve entity, defined as the tag of
   * the Physical Group on that curve whose name is prefixed with
   * "material:" (see cases/README.md)
   *
   * @throw std::runtime_error if no such Physical Group is found
   *
   * @param curve_physical_tags Physical Group tags assigned to the curve
   * @param physical_names Map from Physical Group tag to name
   * @return int
   */
  int GetMaterialID(
      const std::vector<int>& curve_physical_tags,
      const std::unordered_map<int, std::string>& physical_names) const;

  /**
   * @brief Find the source ID for a curve entity, defined as the tag of the
   * Physical Group on the curve whose name is prefixed with "source:" (see
   * cases/README.md)
   *
   * @throw std::runtime_error if no such Physical Group is found
   *
   * @param curve_physical_tags Physical Group tags assigned to the curve
   * @param physical_names Map from Physical Group tag to name
   * @return int
   */
  int GetSourceID(
      const std::vector<int>& curve_physical_tags,
      const std::unordered_map<int, std::string>& physical_names) const;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_MESH_MESH_H_
