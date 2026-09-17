// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#include "problem/sem_problem.h"

#include <gtest/gtest.h>

#include <stdexcept>

#include "mesh/mesh.h"
#include "problem/cg_problem.h"
#include "utils/enums.h"

namespace hummingbird {

TEST(SEMProblemTest, CGFormulationBuildsACGProblem) {
  Mesh mesh;
  mesh.AddNodes({});
  SEMProblem sem_problem(FEFormulation::CG, mesh, /*n_ordinates=*/1);

  ASSERT_NE(sem_problem.get(), nullptr);
  EXPECT_NE(dynamic_cast<CGProblem*>(sem_problem.get()), nullptr);
}

TEST(SEMProblemTest, GetReturnsProblemSizedFromMeshNodeCount) {
  Mesh mesh;
  mesh.AddNodes({Node{}, Node{}, Node{}});
  SEMProblem sem_problem(FEFormulation::CG, mesh, /*n_ordinates=*/1);

  EXPECT_EQ(sem_problem.get()->n_dofs(), mesh.n_nodes());
}

TEST(SEMProblemTest, UnsupportedFormulationThrows) {
  Mesh mesh;
  mesh.AddNodes({});
  const auto unsupported = static_cast<FEFormulation>(-1);

  EXPECT_THROW(SEMProblem(unsupported, mesh, /*n_ordinates=*/1),
              std::runtime_error);
}

}  // namespace hummingbird
