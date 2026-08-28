// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_UTILS_CONSTANTS_H_
#define HUMMINGBIRD_UTILS_CONSTANTS_H_

namespace hummingbird::utils {

/// @brief General tolerance value for double comparisons
constexpr double TOLERANCE = 1e-15;

/// @brief Tolerance value for EXPECT_NEAR in tests
constexpr double EXP_NEAR_TOLERANCE = 1e-14;

}  // namespace hummingbird::utils

#endif  // HUMMINGBIRD_UTILS_CONSTANTS_H_