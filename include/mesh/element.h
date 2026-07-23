#ifndef HUMMINGBIRD_MESH_ELEMENT_H_
#define HUMMINGBIRD_MESH_ELEMENT_H_

#include <memory>
#include <vector>

#include "material.h"
#include "node.h"

namespace hummingbird::mesh {
class Element {
 public:
  Element();

 private:
  const int material_id;
};
}  // namespace hummingbird::mesh

#endif  // HUMMINGBIRD_MESH_ELEMENT_H_
