// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_SIMULATION_H_
#define HUMMINGBIRD_SIMULATION_H_

#include <limits>

namespace hummingbird {
struct Simulation {
  double k_eff = 1.0;
  double scatter_source_l2 = 0;
  double fission_source = 0.0;
  double k_eff_iter_error = std::numeric_limits<double>::infinity();
  double scatter_iter_error = std::numeric_limits<double>::infinity();
};

}  // namespace hummingbird

#endif  // HUMMINGBIRD_SIMULATION_H_
