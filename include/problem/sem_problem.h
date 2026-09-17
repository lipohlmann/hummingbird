// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_PROBLEM_SEM_PROBLEM_H_
#define HUMMINGBIRD_PROBLEM_SEM_PROBLEM_H_

#include "mesh/mesh.h"
#include "problem/cg_problem.h"
#include "problem/problem_base.h"
#include "utils/enums.h"

namespace hummingbird {
class SEMProblem {
 public:
  SEMProblem(const FEFormulation fe_formulation, const Mesh& mesh,
             const size_t n_ordinates);

  ProblemBase* get() const { return sem_problem_.get(); }

 private:
  std::unique_ptr<ProblemBase> sem_problem_;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PROBLEM_SEM_PROBLEM_H_
