// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_PHYSICS_SOURCE_BASE_H_
#define HUMMINGBIRD_PHYSICS_SOURCE_BASE_H_

#include "mesh/node.h"
#include "quadrature/angular/ordinate.h"

namespace hummingbird {
/**
 * @brief Base class for independent sources
 *
 */
class SourceBase {
 public:
  /**
   * @brief Construct a new Source Base object
   *
   */
  SourceBase() = default;

  /**
   * @brief Destroy the Source Base object
   *
   */
  virtual ~SourceBase() = default;

  /**
   * @brief Evaluate the source at a node in a particular direction
   *
   * @param node Node
   * @param ordinate Ordinate (direction)
   * @return Volumetric source in n/cm2/s/str
   */
  virtual double EvaluateAtNode(const Node& node, const Ordinate& ordinate) = 0;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PHYSICS_SOURCE_BASE_H_
