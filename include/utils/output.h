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
 * @brief Print the scalar flux convergence status for the current source
 * iteration
 *
 * @param flux_l2_error L2 norm of the change in scalar flux from the
 * previous iteration
 * @param relative_error Relative L2 error in the scalar flux (flux_l2_error
 * normalized by the current iteration's flux L2 norm)
 * @param iter Iteration number
 */
void print_scatter_status(const double flux_l2_error,
                          const double relative_error, const unsigned int iter);

/**
 * @brief Print a message indicating source iterations have completed
 *
 * @param final_k_eff Final k-eff value
 * @param sim_typ Run mode used for the simulation
 */
void print_scatter_complete(const double final_k_eff, const RunMode sim_typ);

/**
 * @brief Print information about the mesh
 *
 * @param n_nodes Number of nodes in mesh
 * @param n_elements Number of elements in mesh
 * @param dimension Mesh dimension
 */
void print_mesh_info(const size_t n_nodes, const size_t n_elements,
                     const unsigned int dimension);

}  // namespace hummingbird

#endif  // HUMMINGBIRD_OUTPUT_H_
