#ifndef HUMMINGBIRD_MESH_NODE_H_
#define HUMMINGBIRD_MESH_NODE_H_

#include <cstdint>

#include "enums.h"

namespace hummingbird::mesh {
struct Node {
  size_t id;
  double x;
  double y;
  double z;
  Boundary boundary = Boundary::TRANSMISSIVE;
};
}  // namespace hummingbird::mesh

#endif  // HUMMINGBIRD_MESH_NODE_H_
