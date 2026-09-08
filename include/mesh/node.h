// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_MESH_NODE_H_
#define HUMMINGBIRD_MESH_NODE_H_

#include <cstddef>
#include <cstdint>

#include "utils/enums.h"

namespace hummingbird {
struct Node {
  size_t id;
  double x;
  double y;
  double z;
  BC boundary = BC::NONE;

  /// @brief Source ID to access Source Bank
  unsigned int source_id;

  /// @brief Material ID to access Material Bank
  unsigned int material_id;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_MESH_NODE_H_
