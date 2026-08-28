// SPDX-License-Identifier: BSD-3-Clause
//
// Copyright (c) 2026, rcana-project
// Copyright (c) 2026, Liam Pohlmann
//
// Portions of this file are derived from rcana-project and starling.

#ifndef HUMMINGBIRD_UTILS_MISC_H_
#define HUMMINGBIRD_UTILS_MISC_H_

#include "utils/constants.h"

namespace hummingbird::utils {
/**
 * @brief Tests if two double-type numbers are equal using the utils::TOLERANCE
 * value. Taken from https://github.com/starling
 *
 * @param first First number to compare
 * @param second Second number to compare
 * @return true
 * @return false
 */
bool DoubleEqual(const double first, const double second,
                 const double tolerance = TOLERANCE);

/**
 * @brief Compute the relative error between two values
 *
 * @param new_val New value
 * @param old_val Old value
 * @return double
 */
double RelativeError(const double new_val, const double old_val);
}  // namespace hummingbird::utils

#endif  // HUMMINGBIRD_UTILS_MISC_H_
