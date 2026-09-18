#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <memory>
#include <nlohmann/json.hpp>

#include "banks/bc_bank.h"
#include "banks/material_bank.h"
#include "banks/source_bank.h"
#include "mesh/mesh.h"
#include "mesh/node.h"
#include "mesh/segment.h"
#include "problem/cg_problem.h"
#include "quadrature/angular/angular_quadrature_set.h"
#include "quadrature/angular/ordinate.h"
#include "quadrature/gauss_lobatto_legendre.h"
#include "utils/constants.h"

using nlohmann::json;

namespace hummingbird {

namespace {

// Exposes the protected system/forcing/solution vectors so known values can
// be seeded and assembled results inspected directly, without needing a
// full solve cycle. Apply1DBCs itself stays private on CGProblem --
// exercised through the public, non-virtual ProblemBase::ApplyBCs
// dispatcher instead.
class TestableCGProblem : public CGProblem {
 public:
  using CGProblem::CGProblem;
  using ProblemBase::global_forcing_vectors_;
  using ProblemBase::global_system_matrices_;
  using ProblemBase::solution_vectors_;
};

Node MakeNode(size_t id, double x, double y, double z) {
  Node node;
  node.id = id;
  node.x = x;
  node.y = y;
  node.z = z;
  return node;
}

// Builds a MaterialBank containing a single material (looked up by
// GetIDByName("m")) with the given total cross section and zero
// scattering/fission/nu, matching segment_tests.cc's MakeSingleMaterialBank.
MaterialBank MakeSingleMaterialBank(double total_xs) {
  json input = {{"materials",
                 {{"m",
                   {{"scattering_xs", 0.0},
                    {"total_xs", total_xs},
                    {"fission_xs", 0.0},
                    {"nu", 0.0}}}}}};
  return MaterialBank(input);
}

// Builds a MaterialBank with two distinct materials, "mat_a" and "mat_b".
MaterialBank MakeTwoMaterialBank(double total_xs_a, double total_xs_b) {
  json input = {{"materials",
                 {{"mat_a",
                   {{"scattering_xs", 0.0},
                    {"total_xs", total_xs_a},
                    {"fission_xs", 0.0},
                    {"nu", 0.0}}},
                  {"mat_b",
                   {{"scattering_xs", 0.0},
                    {"total_xs", total_xs_b},
                    {"fission_xs", 0.0},
                    {"nu", 0.0}}}}}};
  return MaterialBank(input);
}

// Builds an Ordinate whose x-direction cosine is exactly omega_x, matching
// segment_tests.cc's MakeOrdinateWithXCosine.
Ordinate MakeOrdinateWithXCosine(double omega_x) {
  return Ordinate(std::acos(omega_x), M_PI / 2.0);
}

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

// A SourceBank with no sources beyond the implicit ID-0 "none" entry --
// sufficient for tests that hand-set source_fluxes afterward and never
// resolve a node's source_id away from its default (0).
SourceBank MakeNoSourceBank() { return SourceBank(json{{"sources", json::object()}}); }

// The smallest angular quadrature set AngularQuadratureSet allows (n_polar
// must be >= 2), used only to satisfy InitializeNodeSolutions' GetAbscissa(0)
// call in these single-ordinate (n_ordinates=1) tests.
AngularQuadratureSet MakeMinimalAngularQuad() {
  return AngularQuadratureSet(AngularQuadSet::GL, 1, 2);
}

}  // namespace

// ---------------------------------------------------------------------------
// CGProblem::Apply1DBCs (vacuum), via the public ProblemBase::ApplyBCs
// ---------------------------------------------------------------------------

TEST(CGProblemVacuumBCTest, OutgoingOrdinateAddsToMatrixDiagonal) {
  BCBank bc_bank(MakeSingleVacuumBCJson());
  Mesh mesh = MakeTwoNodeVacuumMesh(bc_bank.GetIDByName("only_bc"));

  // x-direction cosine 0.6 (azimuth=0, polar=asin(0.6)): direction (0.6,0,0.8).
  // At the east node (normal (1,0,0)), Omega.n = 0.6 > 0 -- outgoing.
  Ordinate ordinate(0.0, std::asin(0.6));
  ASSERT_NEAR(ordinate.x(), 0.6, EXP_NEAR_TOLERANCE);

  TestableCGProblem problem(/*n_dofs=*/2, /*n_ordinates=*/1);

  problem.ApplyBCs(mesh, ordinate, bc_bank, 0);

  // Per the SAAF weak form's boundary term (background.tex eq:element-boundary
  // in 2027-ans-mc, moved to the LHS since it involves the unknown psi at
  // this same node/ordinate): matrix(east,east) += Omega.n = 0.6. This is an
  // implicit contribution to the same linear solve, not a forcing term.
  EXPECT_NEAR(problem.global_system_matrices_.at(0)(1, 1), 0.6,
              EXP_NEAR_TOLERANCE);
  // The west node is not outgoing for this ordinate; untouched.
  EXPECT_NEAR(problem.global_system_matrices_.at(0)(0, 0), 0.0,
              EXP_NEAR_TOLERANCE);
  // The forcing vector is never touched by this BC treatment.
  EXPECT_NEAR(problem.global_forcing_vectors_.at(0)(1), 0.0,
              EXP_NEAR_TOLERANCE);
}

TEST(CGProblemVacuumBCTest, IncomingOrdinateAddsNothing) {
  BCBank bc_bank(MakeSingleVacuumBCJson());
  Mesh mesh = MakeTwoNodeVacuumMesh(bc_bank.GetIDByName("only_bc"));

  // Same direction (0.6,0,0.8). At the west node (normal (-1,0,0)),
  // Omega.n = -0.6 < 0 -- incoming, so vacuum means psi = 0 and nothing
  // should be added to the matrix or forcing vector there.
  Ordinate ordinate(0.0, std::asin(0.6));

  TestableCGProblem problem(/*n_dofs=*/2, /*n_ordinates=*/1);

  problem.ApplyBCs(mesh, ordinate, bc_bank, 0);

  EXPECT_NEAR(problem.global_system_matrices_.at(0)(0, 0), 0.0,
              EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(problem.global_forcing_vectors_.at(0)(0), 0.0,
              EXP_NEAR_TOLERANCE);
}

// ---------------------------------------------------------------------------
// Global assembly: AssembleGlobalMatrixData/AssembleGlobalSystem and
// AssembleGlobalForcingData/AssembleGlobalForcing, built directly (no input
// files). Local element matrices/vectors are reused from already-verified
// values in segment_tests.cc; these tests are about assembly (global row/
// col placement, shared-node summation) rather than re-proving local
// element math. Not testing Solve() -- assembly correctness only.
// ---------------------------------------------------------------------------

TEST(CGProblemAssemblyTest,
     TwoUniformElementsThreePointsMatchesHandDerivedGlobalSystem) {
  // Chain of 2 elements, 3 GLL points each, identical length/material:
  // nodes 0,1,2 (element A) and 2,3,4 (element B), node 2 shared. Mesh::
  // Prepare's renumbering (see MeshPrepareTest.PreservesElementConnectivity
  // AfterRenumbering in mesh_tests.cc) gives this exact monotonic
  // left-to-right numbering for a chain built from 3 corner nodes.
  //
  // Reuses the already-verified GLL-3 local stiffness
  // ([[7/6,-4/3,1/6],[-4/3,8/3,-4/3],[1/6,-4/3,7/6]], segment_tests.cc
  // ThreePointQuadraticElementMatchesHandDerivedValues) and mass (weights
  // {1/3,4/3,1/3}, ThreePointMatchesHandDerivedWeights) matrices. With
  // length=4, sigma_t=2, mu_x=1: coeff_K = (1/2)*(2/4) = 1/4, mass diag =
  // sigma_t*length/2*{1/3,4/3,1/3} = {4/3,16/3,4/3}, giving the per-element
  // combined (stiffness+mass) block:
  //   [[13/8, -1/3, 1/24], [-1/3, 6, -1/3], [1/24, -1/3, 13/8]]
  // placed at global {0,1,2} and {2,3,4}; node 2's diagonal sums both
  // elements' (2,2)/(0,0) entries: 13/8 + 13/8 = 13/4.
  //
  // Forcing: source only at the far-right node (global 4, Q=1 elsewhere 0)
  // matches segment_tests.cc's ThreePointSingleBoundarySourceMatchesHand
  // DerivedValues pattern (Q=[0,0,1] in element B's local order), except the
  // derivative term now carries streaming_coeff = mu_x/sigma_t = 1/2 = 0.5
  // (the delta term, i==k, is unaffected): f_B_local =
  // [0.5*1/3*0.5, 0.5*1/3*(-2), 4/2*1/3+0.5*1/3*1.5] = [1/12, -1/3, 11/12].
  // Element A's local source is all zero, so node 2's forcing is untouched by
  // A (0 + 1/12 = 1/12).
  const double length = 4.0;
  const double sigma_t = 2.0;
  MaterialBank material_bank = MakeSingleMaterialBank(sigma_t);
  Ordinate ordinate = MakeOrdinateWithXCosine(1.0);
  GaussLobattoLegendre gll(3);

  Mesh mesh;
  mesh.AddNodes({MakeNode(0, 0.0, 0.0, 0.0), MakeNode(1, length, 0.0, 0.0),
                 MakeNode(2, 2.0 * length, 0.0, 0.0)});
  mesh.AddElement(
      std::make_unique<Segment>(std::array<size_t, 2>{0, 1}, 0, 0, mesh));
  mesh.AddElement(
      std::make_unique<Segment>(std::array<size_t, 2>{1, 2}, 0, 0, mesh));
  mesh.Prepare(gll);
  SourceBank source_bank = MakeNoSourceBank();
  AngularQuadratureSet angular_quad = MakeMinimalAngularQuad();
  mesh.InitializeNodeSolutions(1, *angular_quad.get(), source_bank);
  // Mesh exposes nodes by const reference only; the underlying Mesh (and
  // its nodes) are genuinely non-const here, so mutating a specific
  // already-added node's fields this way is well-defined.
  const_cast<Node&>(mesh.GetNode(4)).source_fluxes.at(0) = 1.0;

  TestableCGProblem problem(/*n_dofs=*/5, /*n_ordinates=*/1);
  auto gmd =
      problem.AssembleGlobalMatrixData(mesh, gll, material_bank, ordinate);
  problem.AssembleGlobalSystem(gmd, 0);
  auto gfd =
      problem.AssembleGlobalForcingData(gll, mesh, material_bank, ordinate, 0);
  problem.AssembleGlobalForcing(gfd, 0);

  const auto& k = problem.global_system_matrices_.at(0);
  ASSERT_EQ(k.n_rows, 5u);
  ASSERT_EQ(k.n_cols, 5u);
  EXPECT_NEAR(k(0, 0), 13.0 / 8.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(0, 1), -1.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(0, 2), 1.0 / 24.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(0, 3), 0.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(0, 4), 0.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(1, 1), 6.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(1, 2), -1.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(2, 2), 13.0 / 4.0, EXP_NEAR_TOLERANCE)
      << "Shared node's diagonal should sum both elements' contributions.";
  EXPECT_NEAR(k(2, 3), -1.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(2, 4), 1.0 / 24.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(3, 3), 6.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(3, 4), -1.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(4, 4), 13.0 / 8.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(1, 0), k(0, 1), EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(2, 0), k(0, 2), EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(2, 1), k(1, 2), EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(3, 2), k(2, 3), EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(4, 2), k(2, 4), EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(4, 3), k(3, 4), EXP_NEAR_TOLERANCE);

  const auto& f = problem.global_forcing_vectors_.at(0);
  ASSERT_EQ(f.n_elem, 5u);
  EXPECT_NEAR(f(0), 0.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(f(1), 0.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(f(2), 1.0 / 12.0, EXP_NEAR_TOLERANCE)
      << "Shared node's forcing should sum both elements' contributions.";
  EXPECT_NEAR(f(3), -1.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(f(4), 11.0 / 12.0, EXP_NEAR_TOLERANCE);
}

TEST(CGProblemAssemblyTest,
     TwoElementsDifferentMaterialsAndLengthsSumSharedNodeCorrectly) {
  // Same 2-element chain topology, but element A (length=2, sigma_t=1) and
  // element B (length=4, sigma_t=2) now differ, so the two local blocks
  // placed into the global matrix are not identical -- a stronger check
  // that assembly places/sums genuinely different per-element contributions
  // rather than coincidentally matching due to symmetry.
  //
  // Element A: coeff_K = (1/1)*(2/2) = 1, mass diag =
  // 1*2/2*{1/3,4/3,1/3} = {1/3,4/3,1/3}:
  //   [[3/2, -4/3, 1/6], [-4/3, 4, -4/3], [1/6, -4/3, 3/2]]
  // Element B: same as the uniform test above:
  //   [[13/8, -1/3, 1/24], [-1/3, 6, -1/3], [1/24, -1/3, 13/8]]
  // Shared node 2: A's right diagonal (3/2) + B's left diagonal (13/8).
  //
  // Forcing: source only at the shared node (global 2, Q=1 elsewhere 0).
  // The derivative term now carries streaming_coeff = mu_x/sigma_t^e (the
  // delta term, i==k, is unaffected by it):
  // For element A (k=2 term, length=2, streaming_coeff=1/1=1, unchanged):
  //   f_A_local = [1/6, -2/3, 5/6].
  // For element B (k=0 term, length=4, streaming_coeff=1/2=0.5):
  //   f_B_local = [2*1/3 + 0.5*1/3*(-1.5), 0.5*1/3*2, 0.5*1/3*(-0.5)]
  //             = [5/12, 1/3, -1/12].
  // Node 2 sums A's right component (5/6) and B's left component (5/12).
  const double length_a = 2.0;
  const double sigma_t_a = 1.0;
  const double length_b = 4.0;
  const double sigma_t_b = 2.0;
  MaterialBank material_bank = MakeTwoMaterialBank(sigma_t_a, sigma_t_b);
  const int mat_a = material_bank.GetIDByName("mat_a");
  const int mat_b = material_bank.GetIDByName("mat_b");
  Ordinate ordinate = MakeOrdinateWithXCosine(1.0);
  GaussLobattoLegendre gll(3);

  Mesh mesh;
  mesh.AddNodes({MakeNode(0, 0.0, 0.0, 0.0), MakeNode(1, length_a, 0.0, 0.0),
                 MakeNode(2, length_a + length_b, 0.0, 0.0)});
  mesh.AddElement(
      std::make_unique<Segment>(std::array<size_t, 2>{0, 1}, mat_a, 0, mesh));
  mesh.AddElement(
      std::make_unique<Segment>(std::array<size_t, 2>{1, 2}, mat_b, 0, mesh));
  mesh.Prepare(gll);
  SourceBank source_bank = MakeNoSourceBank();
  AngularQuadratureSet angular_quad = MakeMinimalAngularQuad();
  mesh.InitializeNodeSolutions(1, *angular_quad.get(), source_bank);
  const_cast<Node&>(mesh.GetNode(2)).source_fluxes.at(0) = 1.0;

  TestableCGProblem problem(/*n_dofs=*/5, /*n_ordinates=*/1);
  auto gmd =
      problem.AssembleGlobalMatrixData(mesh, gll, material_bank, ordinate);
  problem.AssembleGlobalSystem(gmd, 0);
  auto gfd =
      problem.AssembleGlobalForcingData(gll, mesh, material_bank, ordinate, 0);
  problem.AssembleGlobalForcing(gfd, 0);

  const auto& k = problem.global_system_matrices_.at(0);
  EXPECT_NEAR(k(0, 0), 3.0 / 2.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(0, 2), 1.0 / 6.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(1, 1), 4.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(2, 4), 1.0 / 24.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(3, 3), 6.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(4, 4), 13.0 / 8.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(2, 2), 3.0 / 2.0 + 13.0 / 8.0, EXP_NEAR_TOLERANCE)
      << "Shared node's diagonal should sum two different per-element "
         "values.";

  const auto& f = problem.global_forcing_vectors_.at(0);
  EXPECT_NEAR(f(0), 1.0 / 6.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(f(1), -2.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(f(2), 5.0 / 6.0 + 5.0 / 12.0, EXP_NEAR_TOLERANCE)
      << "Shared node's forcing should sum both elements' contributions.";
  EXPECT_NEAR(f(3), 1.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(f(4), -1.0 / 12.0, EXP_NEAR_TOLERANCE);
}

TEST(CGProblemAssemblyTest, ThreeUniformElementsSumBothSharedNodes) {
  // Extends the chain to 3 elements (2 shared nodes: global 2 and 4),
  // checking that placement/summation generalizes past a single shared
  // node. Same per-element block as the first test above (length=4,
  // sigma_t=2, mu_x=1):
  //   [[13/8, -1/3, 1/24], [-1/3, 6, -1/3], [1/24, -1/3, 13/8]]
  const double length = 4.0;
  const double sigma_t = 2.0;
  MaterialBank material_bank = MakeSingleMaterialBank(sigma_t);
  Ordinate ordinate = MakeOrdinateWithXCosine(1.0);
  GaussLobattoLegendre gll(3);

  Mesh mesh;
  mesh.AddNodes({MakeNode(0, 0.0, 0.0, 0.0), MakeNode(1, length, 0.0, 0.0),
                 MakeNode(2, 2.0 * length, 0.0, 0.0),
                 MakeNode(3, 3.0 * length, 0.0, 0.0)});
  mesh.AddElement(
      std::make_unique<Segment>(std::array<size_t, 2>{0, 1}, 0, 0, mesh));
  mesh.AddElement(
      std::make_unique<Segment>(std::array<size_t, 2>{1, 2}, 0, 0, mesh));
  mesh.AddElement(
      std::make_unique<Segment>(std::array<size_t, 2>{2, 3}, 0, 0, mesh));
  mesh.Prepare(gll);

  TestableCGProblem problem(/*n_dofs=*/7, /*n_ordinates=*/1);
  auto gmd =
      problem.AssembleGlobalMatrixData(mesh, gll, material_bank, ordinate);
  problem.AssembleGlobalSystem(gmd, 0);

  const auto& k = problem.global_system_matrices_.at(0);
  ASSERT_EQ(k.n_rows, 7u);
  EXPECT_NEAR(k(0, 0), 13.0 / 8.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(1, 1), 6.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(2, 2), 13.0 / 4.0, EXP_NEAR_TOLERANCE) << "First shared node.";
  EXPECT_NEAR(k(3, 3), 6.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(4, 4), 13.0 / 4.0, EXP_NEAR_TOLERANCE) << "Second shared node.";
  EXPECT_NEAR(k(5, 5), 6.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(6, 6), 13.0 / 8.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(2, 3), -1.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(3, 4), -1.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(2, 4), 1.0 / 24.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(0, 4), 0.0, EXP_NEAR_TOLERANCE)
      << "Non-adjacent elements have no direct coupling.";
  EXPECT_NEAR(k(0, 6), 0.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(k(2, 6), 0.0, EXP_NEAR_TOLERANCE);
}

}  // namespace hummingbird
