// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_PROBLEM_SEM_PROBLEM_H_
#define HUMMINGBIRD_PROBLEM_SEM_PROBLEM_H_

#include "mesh/mesh.h"
#include "problem/cg_problem.h"
#include "problem/problem_base.h"
#include "utils/enums.h"

namespace hummingbird {
/**
 * @brief Wrapper class that owns a ProblemBase of the type corresponding to
 * the given FEFormulation
 *
 */
class SEMProblem {
 public:
  /**
   * @brief Construct a new SEMProblem object, creating the concrete
   * ProblemBase implementation matching fe_formulation
   *
   * @param fe_formulation Finite element formulation to use
   * @param mesh Mesh
   * @param n_ordinates Number of ordinates in angular quadrature
   * @throw std::runtime_error if fe_formulation is not supported
   */
  SEMProblem(const FEFormulation fe_formulation, const Mesh& mesh,
             const size_t n_ordinates);

  /**
   * @brief Get the underlying ProblemBase object
   *
   * @return ProblemBase*
   */
  ProblemBase* get() const { return sem_problem_.get(); }

 private:
  /// @brief Concrete ProblemBase implementation for the chosen formulation
  std::unique_ptr<ProblemBase> sem_problem_;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PROBLEM_SEM_PROBLEM_H_
