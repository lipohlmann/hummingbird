// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_PROBLEM_CG_PROBLEM_H_
#define HUMMINGBIRD_PROBLEM_CG_PROBLEM_H_

#include "problem/problem_base.h"

namespace hummingbird {
/**
 * @brief Class defining a Continuous Galerkin formulation of the SAAF transport
 * equation
 *
 */
class CGProblem : public ProblemBase {
 public:
  /**
   * @brief Construct a new CGProblem object
   *
   * @param n_dofs Number of Degrees of Freedom
   * @param n_ordinates Number of ordinates in angular quadrature
   */
  CGProblem(const size_t n_dofs, const size_t n_ordinates);

  /**
   * @brief See ProblemBase::AssembleGlobalMatrixData. For each element, the
   * local mass and stiffness matrices are summed and scattered into
   * GlobalMatrixData entries keyed by the element's node IDs.
   *
   * @param mesh Mesh
   * @param gll_quad GaussLobattoLegendre quadrature set
   * @param material_bank Material bank
   * @param ordinate Ordinate for equation
   * @return std::vector<GlobalMatrixData>
   */
  std::vector<GlobalMatrixData> AssembleGlobalMatrixData(
      const Mesh& mesh, const GaussLobattoLegendre& gll_quad,
      const MaterialBank& material_bank, const Ordinate& ordinate) override;

  /**
   * @brief See ProblemBase::AssembleGlobalForcingData. For each element, the
   * local forcing vector is scattered into GlobalForcingData entries keyed by
   * the element's node IDs.
   *
   * @param gll_quad GaussLobattoLegendre quadrature set
   * @param mesh Mesh
   * @param material_bank Material bank
   * @param ordinate Ordinate (direction)
   * @param ordinate_index Index in the angular quadrature
   * @return std::vector<GlobalForcingData>
   */
  // this method will need to apply boundary conditions as well if needed
  std::vector<GlobalForcingData> AssembleGlobalForcingData(
      const GaussLobattoLegendre& gll_quad, const Mesh& mesh,
      const MaterialBank& material_bank, const Ordinate& ordinate,
      const size_t ordinate_index) override;

 private:
  /**
   * @brief See ProblemBase::Apply1DBCs. Only vacuum BCs are currently
   * supported.
   *
   * @param mesh Mesh
   * @param ordinate Ordinate (direction)
   * @param bc_bank BCBank object
   * @param ordinate_index Index of the ordinate in the angular quadrature
   * @throw std::runtime_error if a boundary node has a reflective or
   * unsupported BC
   * @todo Need to implement reflective BCs.
   */
  void Apply1DBCs(const Mesh& mesh, const Ordinate& ordinate,
                  const BCBank& bc_bank, const size_t ordinate_index) override;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PROBLEM_CG_PROBLEM_H_
