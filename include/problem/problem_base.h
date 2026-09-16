// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_PROBLEM_PROBLEM_BASE_H_
#define HUMMINGBIRD_PROBLEM_PROBLEM_BASE_H_

#include <armadillo>
#include <vector>

#include "mesh/mesh.h"

namespace hummingbird {
/**
 * @brief Convenience struct for holding the contribution for an entry in a
 * single element to the global matrix.
 *
 */
struct GlobalMatrixData {
  /// @brief Global matrix row ID
  size_t row_id;

  /// @brief Global matrix column ID
  size_t col_id;

  /// @brief Value to be placed in the global matrix. This could be an entry
  /// from a local stiffness or mass matrix, or other depending on the
  /// formulation
  double value;
};

class ProblemBase {
 public:
  /**
   * @brief Assemble the global element system. That is, form \f$ Ax=b\f$.
   *
   * @param global_matrix_data Vector of GlobalMatrixData structs
   */
  void AssembleGlobalSystem(
      const std::vector<GlobalMatrixData>& global_matrix_data);

  /**
   * @brief Solve the linear system using Armadillo's sparse matrix solver.
   * Stores the results in the solution_vector_ member
   *
   */
  void Solve();

 protected:
  /// @brief Global system matrix
  arma::SpMat<double> global_system_matrix_;

  /// @brief Global forcing vector
  arma::Col<double> global_forcing_vector_;

  /// @brief Column vector of the angular flux at the nodes
  arma::Col<double> solution_vector_;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PROBLEM_PROBLEM_BASE_H_
