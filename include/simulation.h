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

  /// @brief L2 norm of the change in scalar flux between the current and
  /// previous source iteration
  double flux_error_l2 = 0;

  /// @brief Fission source
  double fission_source = 0.0;

  /// @brief Relative error in k_eff between the current and previous
  /// iteration
  double k_eff_iter_error = std::numeric_limits<double>::infinity();

  /// @brief Relative L2 error in the scalar flux (flux_error_l2 normalized
  /// by the current iteration's flux L2 norm) between the current and
  /// previous source iteration
  double flux_relative_error = std::numeric_limits<double>::infinity();
};

}  // namespace hummingbird

#endif  // HUMMINGBIRD_SIMULATION_H_
