#ifndef HUMMINGBIRD_MESH_SEGMENT_H_
#define HUMMINGBIRD_MESH_SEGMENT_H_

#include <array>

#include "mesh/element.h"
namespace hummingbird::mesh {
class Segment : public Element {
 private:
  std::array<uint32_t, 2>
      node_indices;  // this might need to change for the quadrature stuff...
};
}  // namespace hummingbird::mesh

#endif  // HUMMINGBIRD_MESH_SEGMENT_H_
