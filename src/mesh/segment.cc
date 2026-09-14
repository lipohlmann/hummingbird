#include "mesh/segment.h"

#include "utils/enums.h"

namespace hummingbird {
Segment::Segment(const std::array<size_t, 2> boundary_node_ids,
                 const int material_id, const int source_id, const Mesh& mesh)
    : boundary_node_ids_(boundary_node_ids),
      length_(ComputeLength(boundary_node_ids, mesh)),
      Element(material_id, source_id) {}

double Segment::ComputeLength(const std::array<size_t, 2> boundary_node_ids,
                              const Mesh& mesh) {
  auto left_node = mesh.GetNode(boundary_node_ids.front());
  auto right_node = mesh.GetNode(boundary_node_ids.back());
  auto distance = DistanceBetweenNodes(left_node, right_node);
  return distance;
}

std::vector<Node> Segment::CreateInteriorNodes(

    const std::vector<Node>& existing_nodes,
    const GaussLobattoLegendre& gll_quadrature) {
  size_t id = existing_nodes.size();

  Node left_node = existing_nodes.at(boundary_node_ids_.at(0));
  Node right_node = existing_nodes.at(boundary_node_ids_.at(1));
  double direction_x = right_node.x - left_node.x;
  double direction_y = right_node.y - left_node.y;
  double direction_z = right_node.z - left_node.z;

  BC interior_boundary = (left_node.boundary == right_node.boundary)
                             ? left_node.boundary
                             : BC::NONE;

  std::vector<Node> nodes(gll_quadrature.n_points() - 2);
  for (auto i = 1; i < gll_quadrature.n_points() - 1; i++) {
    double xi = gll_quadrature.GetAbscissa(i);
    double fraction = (xi + 1.0) / 2.0;
    Node new_node;
    new_node.id = id;
    new_node.x = left_node.x + fraction * direction_x;
    new_node.y = left_node.y + fraction * direction_y;
    new_node.z = left_node.z + fraction * direction_z;
    new_node.boundary = interior_boundary;
    nodes.at(i - 1) = std::move(new_node);
    id++;
  }
  return nodes;
}

arma::SpMat<double> Segment::LocalMassMatrix(
    const GaussLobattoLegendre& gll_quad, const MaterialBank& material_bank) {
  auto sigma_t = material_bank.GetByID(material_id_).total_xs;
  arma::SpMat<double> material_matrix = arma::speye<arma::SpMat<double>>(
      gll_quad.n_points(), gll_quad.n_points());
  for (auto k = 0; k < gll_quad.n_points(); k++)
    material_matrix(k, k) = gll_quad.GetAbscissa(k) * length_ / 2.0 * sigma_t;
  return material_matrix;
}
}  // namespace hummingbird
