// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef STARLING_SIMULATION_H_
#define STARLING_SIMULATION_H_

#include <limits>

namespace starling {
struct Simulation {
  double k_eff = 1.0;
  double scatter_source_l2 = 0;
  double fission_source = 0.0;
  double k_eff_iter_error = std::numeric_limits<double>::infinity();
  double scatter_iter_error = std::numeric_limits<double>::infinity();
};

}  // namespace starling

#endif  // STARLING_SIMULATION_H_
