// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_MESH_NODE_H_
#define HUMMINGBIRD_MESH_NODE_H_

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

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

  /// @brief BC ID to access the boundary condition (BC) bank
  unsigned int bc_id;

  /// @brief Scalar flux on the node
  double scalar_flux = 0;

  /// @brief Angular flux values in order of the SN quadrature set
  std::vector<double> angular_fluxes;

  /// @brief Source flux values in order of the SN quadrature set
  std::vector<double> source_fluxes;
};

/**
 * @brief Compute the Euclidean distance between two nodes
 *
 * @param node_1 First node
 * @param node_2 Second node
 * @return Distance in units of the mesh
 */
inline double DistanceBetweenNodes(const Node& node_1, const Node& node_2) {
  double x_part = node_1.x - node_2.x;
  double y_part = node_1.y - node_2.y;
  double z_part = node_1.z - node_2.z;
  return std::sqrt(x_part * x_part + y_part * y_part + z_part * z_part);
}
}  // namespace hummingbird

#endif  // HUMMINGBIRD_MESH_NODE_H_
