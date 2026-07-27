#ifndef HUMMINGBIRD_MESH_MESH_H_
#define HUMMINGBIRD_MESH_MESH_H_

#include <vector>

#include "mesh/element.h"
#include "mesh/node.h"

namespace hummingbird::mesh {

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
   * @brief Add node to Mesh
   *
   * @param node Node
   */
  void AddNode(const Node node);

  /**
   * @brief Add element to Mesh
   *
   * @param element Element
   */
  void AddElement(const std::unique_ptr<Element> element);

 private:
  /// @brief Nodes in the mesh
  std::vector<Node> nodes_;

  /// @brief Elements in the mesh
  std::vector<std::unique_ptr<Element>> elements_;

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
}  // namespace hummingbird::mesh

#endif  // HUMMINGBIRD_MESH_MESH_H_
