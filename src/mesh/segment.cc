#include "mesh/segment.h"

namespace hummingbird::mesh {
Segment::Segment(const std::array<size_t, 2> boundary_node_ids,
                 const int material_id,
                 const quadrature::GaussLobattoLegendre& gll_quadrature)
    : Element(material_id, gll_quadrature) {}
}  // namespace hummingbird::mesh
