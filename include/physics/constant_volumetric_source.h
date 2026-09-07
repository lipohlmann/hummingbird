// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_PHYSICS_CONSTANT_VOLUMETRIC_SOURCE_H_
#define HUMMINGBIRD_PHYSICS_CONSTANT_VOLUMETRIC_SOURCE_H_

#include "physics/source_base.h"

namespace hummingbird {
/**
 * @brief Defines a constant-value, isotropic volumetric source
 *
 */
class ConstantVolumetricSource : public SourceBase {
 public:
  /**
   * @brief Construct a new Constant Volumetric Source object
   *
   * @param strength Isotropic source strength
   */
  ConstantVolumetricSource(const double strength);

  /**
   * @brief See SourceBase. This overridden function does *not* need any input
   * parameters to be passed.
   *
   */
  double EvaluateAtNode(const Node& node = Node(0, 0, 0, 0),
                        const Ordinate& ordinate = Ordinate(0, 0)) override;

 private:
  /// @brief Isotropic source strength
  const double strength_;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PHYSICS_CONSTANT_VOLUMETRIC_SOURCE_H_
