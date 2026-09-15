#include "mesh/node.h"

namespace hummingbird {
double DistanceBetweenNodes(const Node& node_1, const Node& node_2) {
  double x_part = node_1.x - node_2.x;
  double y_part = node_1.y - node_2.y;
  double z_part = node_1.z - node_2.z;
  return std::sqrt(x_part * x_part + y_part * y_part + z_part * z_part);
}

void UpdateSourceFluxes(Node& node, const MaterialBank& material_bank,
                        const SourceBank& source_bank,
                        const QuadratureBase<Ordinate>& angular_quadrature) {
  for (auto n = 0; n < angular_quadrature.n_points(); n++) {
    double ind_source =
        source_bank.GetByID(node.source_id)
            ->EvaluateAtNode(node, angular_quadrature.GetAbscissa(n));
    double inscattering =
        material_bank.GetByID(node.id).scattering_xs / 2.0 * node.scalar_flux;
    node.source_fluxes[n] = inscattering + ind_source;
  }
}
}  // namespace hummingbird
