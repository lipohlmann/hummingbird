#ifndef HUMMINGBIRD_MESH_SEGMENT_H_
#define HUMMINGBIRD_MESH_SEGMENT_H_

#include <array>

#include "mesh/element.h"
#include "mesh/node.h"
namespace hummingbird::mesh {
class Segment : public Element {
 public:
  /**
   * @brief Construct a new Segment object
   *
   * @param boundary_node_ids Node IDs defining the Segment bounds
   * @param material_id Material ID
   */
  Segment(const std::array<size_t, 2> boundary_node_ids, const int material_id);

  std::vector<Node> CreateInteriorNodes(
      const quadrature::GaussLobattoLegendre& gll_quadrature) override;

 private:
};
}  // namespace hummingbird::mesh

#endif  // HUMMINGBIRD_MESH_SEGMENT_H_
