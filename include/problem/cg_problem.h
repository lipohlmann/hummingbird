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
  std::vector<GlobalMatrixData> AssembleGlobalMatrixData(
      const Mesh& mesh, const GaussLobattoLegendre& gll_quad,
      const MaterialBank& material_bank, const Ordinate& ordinate) override;

  // this method will need to apply boundary conditions as well if needed
  std::vector<GlobalForcingData> AssembleGlobalForcingData(
      const GaussLobattoLegendre& gll_quad, const Mesh& mesh,
      const size_t ordinate_index, const BCBank& bc_bank) override;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PROBLEM_CG_PROBLEM_H_
