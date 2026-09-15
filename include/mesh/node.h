// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_MESH_NODE_H_
#define HUMMINGBIRD_MESH_NODE_H_

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "banks/material_bank.h"
#include "banks/source_bank.h"
#include "quadrature/angular/ordinate.h"
#include "quadrature/quadrature_base.h"
#include "utils/enums.h"

namespace hummingbird {
struct Node {
  /// @brief Unique identifier for the node
  size_t id;

  /// @brief X-coordinate of the node
  double x;

  /// @brief Y-coordinate of the node
  double y;

  /// @brief Z-coordinate of the node
  double z;

  /// @brief Boundary condition type applied at the node
  BC boundary = BC::NONE;

  /// @brief Source ID to access Source Bank
  unsigned int source_id;

  /// @brief Material ID to access Material Bank
  unsigned int material_id;

  /// @brief BC ID to access the boundary condition (BC) bank
  unsigned int bc_id;

  /// @brief Scalar flux on the node
  double scalar_flux;

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
double DistanceBetweenNodes(const Node& node_1, const Node& node_2);

/**
 * @brief Update the scalar flux of a node using the angular quadrature set
 *
 * @param node Node to update
 * @param angular_quadrature Angular quadrature set
 */
void UpdateScalarFlux(Node& node,
                      const QuadratureBase<Ordinate> angular_quadrature);

/**
 * @brief Update the source terms at each angular quadrature point
 *
 * @param node Node to modify
 * @param material_bank Material bank
 * @param source_bank Source bank
 * @param angular_quadrature Angular quadrature set
 */
void UpdateSourceFluxes(Node& node, const MaterialBank& material_bank,
                        const SourceBank& source_bank,
                        const QuadratureBase<Ordinate>& angular_quadrature);
}  // namespace hummingbird

#endif  // HUMMINGBIRD_MESH_NODE_H_
