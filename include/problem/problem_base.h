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

/**
 * @brief Convenience struct for hold the contribution for an entry in a single
 * element to the global forcing vector.
 *
 */
struct GlobalForcingData {
  /// @brief Global forcing vector row ID
  size_t row_id;

  /// @brief Value to be placed in the global vector
  double value;
};

class ProblemBase {
 public:
  /**
   * @brief Assemble the global element system. That is, form the matrix \f$A\f$
   * in \f$ Ax=b\f$. This only should be called once per ordinate per simulation
   *
   * @param global_matrix_data Vector of GlobalMatrixData structs
   */
  void AssembleGlobalSystem(
      const std::vector<GlobalMatrixData>& global_matrix_data,
      const size_t ordinate_index);

  /**
   * @brief Solve the linear system using Armadillo's sparse matrix solver.
   * Stores the results in the solution_vectors_ member at the ordinate_index
   * index
   *
   */
  void Solve(const size_t ordinate_index);

  /**
   * @brief Assemble the vector of GlobalMatrixData structs that will be used in
   * AssembleGlobalSystem to create the global system matrix.
   *
   * @param mesh Mesh
   * @param gll_quad GaussLobattoLegendre quadrature set
   * @param material_bank Material bank
   * @param ordinate Ordinate for equation
   * @return std::vector<GlobalMatrixData>
   */
  virtual std::vector<GlobalMatrixData> AssembleGlobalMatrixData(
      const Mesh& mesh, const GaussLobattoLegendre& gll_quad,
      const MaterialBank& material_bank, const Ordinate& ordinate) = 0;

 protected:
  /// @brief Global system matrices in order of ordinates in angular quadrature
  /// set
  std::vector<arma::SpMat<double>> global_system_matrices_;

  /// @brief Global forcing vectors in order of ordinates in angular quadrature
  /// set
  std::vector<arma::Col<double>> global_forcing_vectors_;

  /// @brief Column vector of the angular flux at the nodes in order of
  /// ordinates in angular quadrature set
  std::vector<arma::Col<double>> solution_vectors_;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PROBLEM_PROBLEM_BASE_H_
