#ifndef HUMMINGBIRD_NODE_H_
#define HUMMINGBIRD_NODE_H_

#include <cstdint>

namespace hummingbird {
struct Node {
  uint32_t id;
  double x;
  double y;
  double z;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_NODE_H_
