#ifndef HUMMINGBIRD_ELEMENT_H
#define HUMMINGBIRD_ELEMENT_H

#include <memory>
#include <vector>

#include "material.h"
#include "node.h"

namespace hummingbird {
class Element {
 public:
  Element();

 private:
  const int material_id;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_ELEMENT_H
