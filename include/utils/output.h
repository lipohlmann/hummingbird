// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_OUTPUT_H_
#define HUMMINGBIRD_OUTPUT_H_

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <string>

#include "utils/enums.h"

namespace hummingbird {

/**
 * @brief Print the program banner, description, and license information
 *
 */
void print_header();

/**
 * @brief Print the paths of the input and mesh files being used
 *
 * @param input Path to the input file
 * @param mesh Path to the mesh file
 */
void print_input_files(const std::string input, const std::string mesh);

/**
 * @brief Print the column headers for the k-eff/scatter iteration status
 * table
 *
 */
void print_columns();

/**
 * @brief Print the k-eff status for the current iteration
 *
 * @param k_eff Current k-eff value
 * @param error Relative error in k-eff from the previous iteration
 * @param iter Iteration number
 */
void print_k_status(const double k_eff, const double error,
                    const unsigned int iter);

/**
 * @brief Print the scattering source status for the current iteration
 *
 * @param scatter Current scattering source L2 norm
 * @param error Relative error in the scattering source from the previous
 * iteration
 * @param iter Iteration number
 */
void print_scatter_status(const double scatter, const double error,
                          const unsigned int iter);

/**
 * @brief Print a message indicating source iterations have completed
 *
 * @param final_k_eff Final k-eff value
 * @param sim_typ Run mode used for the simulation
 */
void print_scatter_complete(const double final_k_eff, const RunMode sim_typ);

}  // namespace hummingbird

#endif  // HUMMINGBIRD_OUTPUT_H_
