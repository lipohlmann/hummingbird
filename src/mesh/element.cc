#include "mesh/element.h"

#include <format>
#include <stdexcept>

#include "mesh/mesh.h"

namespace hummingbird {
Element::Element(const int material_id, const int source_id)
    : material_id_(material_id), source_id_(source_id) {}

void Element::SetNewNodeID(const size_t prev_id, const size_t new_id) {
  for (auto i = 0; i < node_ids_.size(); i++) {
    if (node_ids_.at(i) == prev_id) {
      node_ids_.at(i) = new_id;
      return;
    }
  }
  throw std::runtime_error(
      std::format("ID {} not found in current element.", prev_id));
}

void Element::SetNodeSourceFluxes(
    Mesh& mesh, const MaterialBank& material_bank,
    const SourceBank& source_bank,
    const QuadratureBase<Ordinate>& angular_quadrature) {
  for (const auto node_id : node_ids_) {
    const Node& node = mesh.GetNode(node_id);
    auto source_fluxes =
        ComputeNodeSourceFluxes(node, material_id_, source_id_, material_bank,
                                source_bank, angular_quadrature);
    mesh.SetNodeSourceFluxes(node_id, source_fluxes);
  }
}
}  // namespace hummingbird
