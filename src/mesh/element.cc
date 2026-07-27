#include "mesh/element.h"

#include <format>
#include <stdexcept>

namespace hummingbird::mesh {
Element::Element(const int material_id,
                 const quadrature::GaussLobattoLegendre& gll_quadrature)
    : material_id_(material_id), gll_quadrature_(gll_quadrature) {
  node_ids_.reserve(gll_quadrature.n_points());
}

void Element::SetNewNodeID(const size_t prev_id, const size_t new_id) {
  for (auto i = 0; i < node_ids_.size(); i++)
    if (node_ids_.at(i) == prev_id) node_ids_.at(i) = new_id;
  throw std::runtime_error(
      std::format("ID {} not found in current element.", prev_id));
}
}  // namespace hummingbird::mesh
