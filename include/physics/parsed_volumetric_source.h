// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_PHYSICS_PARSED_VOLUMETRIC_SOURCE_H_
#define HUMMINGBIRD_PHYSICS_PARSED_VOLUMETRIC_SOURCE_H_

#include <string>

#include "mesh/node.h"
#include "physics/source_base.h"
#include "quadrature/angular/ordinate.h"

namespace hummingbird {
/**
 * @brief Defines a volumetric source defined by a parsed function passed by the
 * user
 *
 */
class ParsedVolumetricSource : public SourceBase {
 public:
  /**
   * @brief Construct a new Parsed Volumetric Source object
   *
   * @param expression Function expression to be evaluated
   */
  ParsedVolumetricSource(const std::string expression);

  /**
   * @brief Evaluate the parsed expression_ at a node in a particular
   * direction. The expression may reference the variables x, y, z (the
   * node's coordinates) and mu (the direction's x-direction cosine)
   *
   * @param node Node
   * @param ordinate Ordinate (direction)
   * @return Volumetric source in n/cm2/s/str
   * @throw std::runtime_error if the expression fails to parse
   * @todo Add support for 2- and 3D angular functions (currently only supports
   * mu-functions).
   */
  double EvaluateAtNode(const Node& node, const Ordinate& ordinate) override;

 private:
  /// @brief Function expression to be evaluated
  const std::string expression_;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PHYSICS_PARSED_VOLUMETRIC_SOURCE_H_
