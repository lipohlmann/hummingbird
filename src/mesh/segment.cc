#include "mesh/segment.h"

namespace hummingbird::mesh {
Segment::Segment(const std::array<size_t, 2> boundary_node_ids,
                 const int material_id)
    : Element(material_id) {}

std::vector<Node> Segment::CreateInteriorNodes(
    const quadrature::GaussLobattoLegendre& gll_quadrature) {
  node_ids_.reserve(gll_quadrature.n_points());
}
}  // namespace hummingbird::mesh
