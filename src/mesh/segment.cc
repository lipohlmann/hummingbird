#include "mesh/segment.h"

#include "enums.h"

namespace hummingbird {
Segment::Segment(const std::array<size_t, 2> boundary_node_ids,
                 const int material_id)
    : boundary_node_ids_(boundary_node_ids), Element(material_id) {}

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
}  // namespace hummingbird
