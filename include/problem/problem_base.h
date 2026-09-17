// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_PROBLEM_PROBLEM_BASE_H_
#define HUMMINGBIRD_PROBLEM_PROBLEM_BASE_H_

#include <armadillo>
#include <vector>

#include "banks/bc_bank.h"
#include "banks/material_bank.h"
#include "mesh/mesh.h"
#include "quadrature/quadrature_base.h"

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

/**
 * @brief Base class defining a formulation of the transport equation to be
 * solved
 *
 */
class ProblemBase {
 public:
  /**
   * @brief Construct a new Problem Base object
   *
   * @param n_dofs Number of Degrees of Freedom
   * @param n_ordinates Number of ordinates in angular mesh
   */
  ProblemBase(const size_t n_dofs, const size_t n_ordinates);

  /**
   * @brief Destroy the Problem Base object
   *
   */
  virtual ~ProblemBase() = default;

  /**
   * @brief Assemble the global element system. That is, form the matrix \f$A\f$
   * in \f$ Ax=b\f$. This only should be called once per ordinate per
   * simulation. The global matrix is sized from the n_dofs supplied at
   * construction, not from global_matrix_data; row_id/col_id pairs are
   * expected to repeat (shared nodes get summed contributions from multiple
   * elements) and are summed via Armadillo's batch SpMat constructor's
   * add_values=true overload (the default overload errors on duplicate
   * locations instead of summing them).
   *
   * @param global_matrix_data Vector of GlobalMatrixData structs
   * @param ordinate_index Ordinate index within the angular quadrature set
   * (i.e. the direction number)
   */
  void AssembleGlobalSystem(
      const std::vector<GlobalMatrixData>& global_matrix_data,
      const size_t ordinate_index);

  /**
   * @brief Assemble the global forcing vector for a given ordinate index. The
   * global vector is sized from the n_dofs supplied at construction, not from
   * global_forcing_data; row_id values are expected to repeat (shared nodes
   * get summed contributions from multiple elements).
   *
   * @param global_forcing_data Vector of GlobalForcingData structs
   * @param ordinate_index Ordinate index
   */
  void AssembleGlobalForcing(
      const std::vector<GlobalForcingData>& global_forcing_data,
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

  /**
   * @brief Assemble the vector of GlobalForcingData structs that will be using
   * in AssembleGlobalForcing to create the global forcing vector
   *
   * @param gll_quad GaussLobattoLegendre quadrature set
   * @param mesh Mesh
   * @param ordinate_index Index in the angular quadrature
   * @return std::vector<GlobalForcingData>
   */
  virtual std::vector<GlobalForcingData> AssembleGlobalForcingData(
      const GaussLobattoLegendre& gll_quad, const Mesh& mesh,
      const size_t ordinate_index) = 0;

  /**
   * @brief Apply boundary conditions by adding in values needed at boundary
   * nodes. Internally, a switch calls the correct methods depending on the mesh
   * dimension.
   *
   * @param mesh Mesh
   * @param ordinate Ordinate (direction)
   * @param bc_bank BCBank object
   * @param ordinate_index Index of the ordinate in the angular quadrature
   */
  void ApplyBCs(const Mesh& mesh, const Ordinate& ordinate,
                const BCBank& bc_bank, const size_t ordinate_index);

  /**
   * @brief Get the solution vectors
   *
   * @return const std::vector<arma::Col<double>>&
   */
  const std::vector<arma::Col<double>>& solution_vectors() const {
    return solution_vectors_;
  }

  /**
   * @brief Get the number of degrees of freedom
   *
   * @return size_t
   */
  size_t n_dofs() const { return n_dofs_; }

 protected:
  /// @brief Number of degrees of freedom in simulation
  const size_t n_dofs_;

  /// @brief Global system matrices in order of ordinates in angular
  /// quadrature set
  std::vector<arma::SpMat<double>> global_system_matrices_;

  /// @brief Global forcing vectors in order of ordinates in angular quadrature
  /// set
  std::vector<arma::Col<double>> global_forcing_vectors_;

  /// @brief Column vector of the angular flux at the nodes in order of
  /// ordinates in angular quadrature set
  std::vector<arma::Col<double>> solution_vectors_;

  /**
   * @brief Validate a vector of GlobalMatrixData before it is used to
   * assemble the global system matrix. row_id/col_id pairs are expected to
   * repeat (shared nodes get summed contributions from multiple elements),
   * so unique (row_id, col_id) pairs are not required. Instead, this checks
   * that the set of distinct row/col IDs referenced is contiguous starting
   * from 0, i.e. no node ID is missing from the data.
   *
   * @param gmd Vector of GlobalMatrixData structs
   */
  void CheckGlobalMatrixData(const std::vector<GlobalMatrixData>& gmd);

  /**
   * @brief Validate a vector of GlobalForcingData before it is used to
   * assemble the global forcing vector. row_id values are expected to repeat
   * (shared nodes get summed contributions from multiple elements), so
   * unique row_ids are not required. Instead, this checks that the set of
   * distinct row IDs referenced is contiguous starting from 0, i.e. no DOF
   * is missing from the data.
   *
   * @param gfd Vector of GlobalForcingData structs
   */
  void CheckGlobalForcingData(const std::vector<GlobalForcingData>& gfd);

  /**
   * @brief Apply boundary conditions to a 1D problem.
   *
   * @param mesh Mesh
   * @param Ordinate Ordinate (direction)
   * @param bc_bank BCBank object
   * @param ordinate_index Index of the ordinate in the angular quadrature
   */
  virtual void Apply1DBCs(const Mesh& mesh, const Ordinate& Ordinate,
                          const BCBank& bc_bank,
                          const size_t ordinate_index) = 0;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PROBLEM_PROBLEM_BASE_H_
