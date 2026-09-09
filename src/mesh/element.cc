#include "mesh/element.h"

#include <format>
#include <stdexcept>

namespace hummingbird {
Element::Element(const int material_id, const int source_id)
    : material_id_(material_id), source_id_(source_id) {}

void Element::SetNewNodeID(const size_t prev_id, const size_t new_id) {
  for (auto i = 0; i < node_ids_.size(); i++)
    if (node_ids_.at(i) == prev_id) node_ids_.at(i) = new_id;
  throw std::runtime_error(
      std::format("ID {} not found in current element.", prev_id));
}
}  // namespace hummingbird
