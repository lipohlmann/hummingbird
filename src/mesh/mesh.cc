#include "mesh/mesh.h"

#include <algorithm>

namespace hummingbird::mesh {
void Mesh::RenumberNodes() {
  std::sort(nodes_.begin(), nodes_.end(),
            [](const Node& first, const Node& second) {
              return first.id < second.id;
            });

  uint32_t new_node_id = 0;
  for (auto& element : elements_) {
    auto node_ids = element->node_ids();
    for (auto id : node_ids) {
      // Elements will share nodes. Don't renumber if we've already set a new ID
      if (id < new_node_id) continue;
      nodes_.at(id).id = new_node_id;
      new_node_id++;
    }
  }
}
}  // namespace hummingbird::mesh
