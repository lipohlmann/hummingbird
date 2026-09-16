#include "problem/cg_problem.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <memory>
#include <nlohmann/json.hpp>

#include "banks/bc_bank.h"
#include "mesh/mesh.h"
#include "mesh/node.h"
#include "mesh/segment.h"
#include "quadrature/angular/ordinate.h"
#include "utils/constants.h"

using nlohmann::json;

namespace hummingbird {

namespace {

// Exposes the protected forcing/solution vectors so a known lagged solution
// value can be seeded and the resulting forcing vector inspected, without
// needing to drive a full assemble/solve cycle. Apply1DBCs itself stays
// private on CGProblem -- exercised through the public, non-virtual
// ProblemBase::ApplyBCs dispatcher instead.
class TestableCGProblem : public CGProblem {
 public:
  using CGProblem::CGProblem;
  using ProblemBase::global_forcing_vectors_;
  using ProblemBase::solution_vectors_;
};

// A minimal 1D, 2-node mesh (one Segment) with both endpoints tagged VACUUM
// and their outward normals set by hand, bypassing gmsh/ResolveIDs so the
// BC ids and geometry under test are explicit and self-contained.
Mesh MakeTwoNodeVacuumMesh(unsigned int vacuum_bc_id) {
  Node west;
  west.id = 0;
  west.x = 0.0;
  west.y = 0.0;
  west.z = 0.0;
  west.bc_id = vacuum_bc_id;
  west.outward_normal = {-1.0, 0.0, 0.0};

  Node east;
  east.id = 1;
  east.x = 1.0;
  east.y = 0.0;
  east.z = 0.0;
  east.bc_id = vacuum_bc_id;
  east.outward_normal = {1.0, 0.0, 0.0};

  Mesh mesh;
  mesh.AddNodes({west, east});
  mesh.AddElement(
      std::make_unique<Segment>(std::array<size_t, 2>{0, 1}, 0, 0, mesh));
  mesh.FindBoundaryNodes();
  return mesh;
}

json MakeSingleVacuumBCJson() {
  return json{{"boundary_conditions", {{"only_bc", {{"type", "vacuum"}}}}}};
}

}  // namespace

// ---------------------------------------------------------------------------
// CGProblem::Apply1DBCs (vacuum), via the public ProblemBase::ApplyBCs
// ---------------------------------------------------------------------------

TEST(CGProblemVacuumBCTest, OutgoingOrdinateAddsNegatedLaggedSolution) {
  BCBank bc_bank(MakeSingleVacuumBCJson());
  Mesh mesh = MakeTwoNodeVacuumMesh(bc_bank.GetIDByName("only_bc"));

  // x-direction cosine 0.6 (azimuth=0, polar=asin(0.6)): direction (0.6,0,0.8).
  // At the east node (normal (1,0,0)), Omega.n = 0.6 > 0 -- outgoing.
  Ordinate ordinate(0.0, std::asin(0.6));
  ASSERT_NEAR(ordinate.x(), 0.6, EXP_NEAR_TOLERANCE);

  TestableCGProblem problem(/*n_dofs=*/2, /*n_ordinates=*/1);
  problem.solution_vectors_.at(0)(1) = 5.0;  // lagged value at the east node

  problem.ApplyBCs(mesh, ordinate, bc_bank, 0);

  // Expected: forcing(east) += -(Omega.n) * solution(east) = -0.6 * 5.0
  EXPECT_NEAR(problem.global_forcing_vectors_.at(0)(1), -3.0,
             EXP_NEAR_TOLERANCE);
  // The west node is not outgoing for this ordinate; untouched.
  EXPECT_NEAR(problem.global_forcing_vectors_.at(0)(0), 0.0,
             EXP_NEAR_TOLERANCE);
}

TEST(CGProblemVacuumBCTest, IncomingOrdinateAddsNothing) {
  BCBank bc_bank(MakeSingleVacuumBCJson());
  Mesh mesh = MakeTwoNodeVacuumMesh(bc_bank.GetIDByName("only_bc"));

  // Same direction (0.6,0,0.8). At the west node (normal (-1,0,0)),
  // Omega.n = -0.6 < 0 -- incoming, so vacuum means psi = 0 and nothing
  // should be added to the forcing vector there.
  Ordinate ordinate(0.0, std::asin(0.6));

  TestableCGProblem problem(/*n_dofs=*/2, /*n_ordinates=*/1);
  problem.solution_vectors_.at(0)(0) = 5.0;  // lagged value at the west node

  problem.ApplyBCs(mesh, ordinate, bc_bank, 0);

  EXPECT_NEAR(problem.global_forcing_vectors_.at(0)(0), 0.0,
             EXP_NEAR_TOLERANCE);
}

}  // namespace hummingbird
