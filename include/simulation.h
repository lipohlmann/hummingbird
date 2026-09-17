// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_SIMULATION_H_
#define HUMMINGBIRD_SIMULATION_H_

#include <limits>

namespace hummingbird {

/**
 * @brief Tracks the running state of a simulation across source iterations
 *
 */
struct Simulation {
  /// @brief Effective multiplication factor
  double k_eff = 1.0;

  /// @brief L2 norm of the scattering source
  double scatter_source_l2 = 0;

  /// @brief Fission source
  double fission_source = 0.0;

  /// @brief Relative error in k_eff between the current and previous
  /// iteration
  double k_eff_iter_error = std::numeric_limits<double>::infinity();

  /// @brief Relative error in the scattering source L2 norm between the
  /// current and previous iteration
  double scatter_iter_error = std::numeric_limits<double>::infinity();
};

}  // namespace hummingbird

#endif  // HUMMINGBIRD_SIMULATION_H_
