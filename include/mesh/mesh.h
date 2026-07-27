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
  Mesh();

 private:
  std::vector<Node> nodes_;
  std::vector<std::unique_ptr<Element>> elements_;

  /**
   * @brief Renumber nodes in mesh to keep node IDs near each other in a single
   * element
   *
   */
  void RenumberNodes();
};
}  // namespace hummingbird::mesh

#endif  // HUMMINGBIRD_MESH_MESH_H_
