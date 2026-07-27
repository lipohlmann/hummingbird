#include "mesh/element.h"

namespace hummingbird::mesh {
Element::Element(const int material_id,
                 const quadrature::GaussLobattoLegendre& gll_quadrature)
    : material_id_(material_id), gll_quadrature_(gll_quadrature) {
  nodes_.reserve(gll_quadrature.n_points());
}
}  // namespace hummingbird::mesh
