// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_MESH_MESH_H_
#define HUMMINGBIRD_MESH_MESH_H_

#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "banks/bc_bank.h"
#include "banks/material_bank.h"
#include "banks/source_bank.h"
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
   * @brief Get the elements in the mesh
   *
   * @return const std::vector<std::unique_ptr<Element>>&
   */
  const std::vector<std::unique_ptr<Element>>& elements() const {
    return elements_;
  }

  /**
   * @brief Get Node by ID
   *
   * @param id NodeID
   * @return const Node&
   */
  const Node& GetNode(const size_t id) const { return nodes_.at(id); }

  /**
   * @brief Get the number of Elements in the mesh
   *
   * @return size_t
   */
  size_t n_elements() const { return elements_.size(); }

  /**
   * @brief Get the number of nodes in the mesh
   *
   * @return size_t
   */
  size_t n_nodes() const { return nodes_.size(); }

  /**
   * @brief Get the spatial dimension of the mesh, derived from the elements
   * it contains (e.g. 1 for a mesh of Segments). 0 if no elements have been
   * added yet.
   *
   * @return unsigned int
   */
  unsigned int dimension() const { return dimension_; }

  /**
   * @brief Get an Element by index
   *
   * @param index Index of the element
   * @return const Element&
   */
  const Element& GetElement(const size_t index) const {
    return *elements_.at(index);
  }

  /**
   * @brief Initialize all node solution values to 0's with the correct vector
   * lengths
   *
   * @param n_ordinates Number of ordinates
   */
  void InitializeNodeSolutions(const size_t n_ordinates);

  /**
   * @brief Prepare the mesh for running a simulation by generating interior
   * nodes using the GLL quadrature set, renumbering nodes, then checking these
   * nodes.
   *
   * @param gll_quadrature GaussLegendreLobatto quadrature set
   */
  void Prepare(const GaussLobattoLegendre& gll_quadrature);

  /**
   * @brief Resolve the raw gmsh Physical Group tags currently stored in
   * material_id()/source_id()/bc_id (set from the .msh file alone, with no
   * bank access) into the real IDs assigned by the given banks, by looking
   * up the name each tag's Physical Group carries (e.g. "mms_material" from
   * "material:mms_material") via BankBase::GetIDByName. Nodes with
   * bc_id == 0 (never tagged by a point element) are left alone.
   *
   * @throw std::runtime_error if a material name from the mesh has no match
   * in material_bank
   *
   * @param material_bank Material bank
   * @param source_bank Source bank
   * @param bc_bank BC bank
   */
  void ResolveIDs(const MaterialBank& material_bank,
                  const SourceBank& source_bank, const BCBank& bc_bank);

  /**
   * @brief Populates the boundary_node_ids_ member. Must be called after
   * ResolveIDs!
   *
   */
  void FindBoundaryNodes();

 private:
  /// @brief Nodes in the mesh
  std::vector<Node> nodes_;

  /// @brief Nodes in the mesh that exist on the mesh boundaries (that is, they
  /// have a boundary condition assigned)
  std::vector<size_t> boundary_node_ids_;

  /// @brief Elements in the mesh
  std::vector<std::unique_ptr<Element>> elements_;

  /// @brief Spatial dimension of the mesh's elements, set by the first call
  /// to AddElement. 0 if no elements have been added yet.
  unsigned int dimension_ = 0;

  /// @brief Map from Physical Group tag to its name (e.g. "material:foo"),
  /// persisted from GmshReadState::physical_names for use by ResolveIDs
  /// after reading is done.
  std::unordered_map<int, std::string> physical_names_;

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
   * for 2-node line elements. Point elements are consumed and used to tag
   * the boundary/bc_id of the Node they reference, since they mark boundary
   * entities in gmsh.
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

  /**
   * @brief Find the BC ID for a point entity, defined as the tag of the
   * Physical Group on that point whose name is prefixed with "bc:" (see
   * cases/README.md)
   *
   * @throw std::runtime_error if no such Physical Group is found
   *
   * @param point_physical_tags Physical Group tags assigned to the point
   * @param physical_names Map from Physical Group tag to name
   * @return int
   */
  int GetBCID(const std::vector<int>& point_physical_tags,
              const std::unordered_map<int, std::string>& physical_names) const;

  /**
   * @brief Extract the name after the ":" in a Physical Group name (e.g.
   * "mms_material" from "material:mms_material")
   *
   * @param physical_name Physical Group name
   * @return std::string
   */
  std::string ExtractName(const std::string& physical_name) const;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_MESH_MESH_H_
