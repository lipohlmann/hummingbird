#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "banks/material_bank.h"
#include "mesh/mesh.h"
#include "mesh/node.h"
#include "mesh/segment.h"
#include "physics/material.h"
#include "quadrature/angular/ordinate.h"
#include "quadrature/gauss_lobatto_legendre.h"
#include "utils/constants.h"

using nlohmann::json;

namespace hummingbird {

namespace {

Node MakeNode(size_t id, double x, double y, double z, unsigned int bc_id = 0) {
  Node node;
  node.id = id;
  node.x = x;
  node.y = y;
  node.z = z;
  node.bc_id = bc_id;
  return node;
}

// Builds a MaterialBank containing a single material (ID 0) with the given
// total cross section and zero scattering/fission/nu, so LocalStiffnessMatrix
// and LocalMassMatrix tests can isolate the effect of total_xs.
MaterialBank MakeSingleMaterialBank(double total_xs) {
  json input = {{"materials",
                 {{"test_material",
                   {{"scattering_xs", 0.0},
                    {"total_xs", total_xs},
                    {"fission_xs", 0.0},
                    {"nu", 0.0}}}}}};
  return MaterialBank(input);
}

// Builds an Ordinate whose x-direction cosine is exactly omega_x, by fixing
// the polar angle at pi/2 (so sin(polar) == 1) and solving azimuth =
// acos(omega_x).
Ordinate MakeOrdinateWithXCosine(double omega_x) {
  return Ordinate(std::acos(omega_x), M_PI / 2.0);
}

}  // namespace

// ---------------------------------------------------------------------------
// Interior node count
// ---------------------------------------------------------------------------

class SegmentInteriorCountTest : public ::testing::TestWithParam<size_t> {
 protected:
  Mesh mesh;
};

TEST_P(SegmentInteriorCountTest, ProducesNPointsMinusTwoInteriorNodes) {
  const size_t n_points = GetParam();
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 10.0, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 1, 0, mesh);
  GaussLobattoLegendre gll(n_points);

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  EXPECT_EQ(interior.size(), n_points - 2);
}

INSTANTIATE_TEST_SUITE_P(VariousOrders, SegmentInteriorCountTest,
                         ::testing::Values(2, 3, 4, 5, 6, 8));

class SegmentTest : public testing::Test {
 protected:
  Mesh mesh;
};

TEST_F(SegmentTest, TwoPointQuadratureProducesNoInteriorNodes) {
  // A GLL rule with only the two endpoints has no interior points at all.
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 10.0, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 1, 0, mesh);
  GaussLobattoLegendre gll(2);

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  EXPECT_TRUE(interior.empty());
}

TEST_F(SegmentTest, DimensionIsOne) {
  mesh.AddNodes({MakeNode(0, 0.0, 0.0, 0.0), MakeNode(1, 10.0, 0.0, 0.0)});
  Segment segment({0, 1}, 1, 0, mesh);
  EXPECT_EQ(segment.dimension(), 1u);
}

// ---------------------------------------------------------------------------
// Node ID assignment
// ---------------------------------------------------------------------------

TEST_F(SegmentTest, InteriorNodeIdsContinueFromExistingNodeCount) {
  // Three nodes already exist, so new IDs should start at 3 and increment.
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 10.0, 0.0, 0.0),
                             MakeNode(2, 20.0, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 1, 0, mesh);
  GaussLobattoLegendre gll(5);  // 3 interior points

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_EQ(interior.size(), 3u);
  for (size_t k = 0; k < interior.size(); ++k) {
    EXPECT_EQ(interior.at(k).id, nodes.size() + k)
        << "Interior node " << k << " has an unexpected ID.";
  }
}

TEST_F(SegmentTest, UsesNodesAtSpecifiedBCIndicesRegardlessOfPosition) {
  // boundary_node_ids_ should be used as indices into existing_nodes, not
  // assumed to be the first two entries.
  std::vector<Node> nodes(6);
  nodes.at(2) = MakeNode(2, 0.0, 0.0, 0.0);
  nodes.at(5) = MakeNode(5, 10.0, 0.0, 0.0);
  mesh.AddNodes(nodes);
  Segment segment({2, 5}, 1, 0, mesh);
  GaussLobattoLegendre gll(3);  // 1 interior point

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_EQ(interior.size(), 1u);
  EXPECT_EQ(interior.at(0).id, nodes.size());
}

TEST_F(SegmentTest, CreateInteriorNodesPopulatesNodeIds) {
  // node_ids() should end up as [left, interior IDs in order, right],
  // regardless of the geometric mapping.
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 10.0, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 1, 0, mesh);
  GaussLobattoLegendre gll(5);  // 3 interior points -> IDs 2, 3, 4

  EXPECT_TRUE(segment.node_ids().empty());
  auto interior = segment.CreateInteriorNodes(nodes, gll);

  ASSERT_EQ(interior.size(), 3u);
  EXPECT_EQ(segment.node_ids(), (std::vector<size_t>{0, 2, 3, 4, 1}));
}

TEST_F(SegmentTest, CreateInteriorNodesNodeIdsFollowBoundaryOrder) {
  // Reversing boundary_node_ids_ should reverse which endpoint's ID comes
  // first in node_ids(), matching DirectionFollowsBCNodeIdOrder above.
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 10.0, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({1, 0}, 1, 0, mesh);
  GaussLobattoLegendre gll(3);  // 1 interior point -> ID 2

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_EQ(interior.size(), 1u);
  EXPECT_EQ(segment.node_ids(), (std::vector<size_t>{1, 2, 0}));
}

// ---------------------------------------------------------------------------
// Coordinate mapping
//
// Per Segment::CreateInteriorNodes's documented mapping,
//   r(xi) = r1 + ((xi + 1) / 2) * (r2 - r1),  xi in [-1, 1]
// where r1 is the node at boundary_node_ids_[0] and r2 is the node at
// boundary_node_ids_[1].
// ---------------------------------------------------------------------------

TEST_F(SegmentTest, InteriorNodesMatchDocumentedReferenceToPhysicalMapping) {
  const Node left = MakeNode(0, 1.0, 2.0, 3.0);
  const Node right = MakeNode(1, 4.0, 6.0, 3.0);  // non-axis-aligned segment
  std::vector<Node> nodes = {left, right};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 1, 0, mesh);
  GaussLobattoLegendre gll(6);  // 4 interior points

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_EQ(interior.size(), gll.n_points() - 2);

  for (size_t k = 0; k < interior.size(); ++k) {
    const size_t gll_index = k + 1;  // matches the implementation's loop
    const double xi = gll.GetAbscissa(gll_index);
    const double fraction = (xi + 1.0) / 2.0;

    const double expected_x = left.x + fraction * (right.x - left.x);
    const double expected_y = left.y + fraction * (right.y - left.y);
    const double expected_z = left.z + fraction * (right.z - left.z);

    EXPECT_NEAR(interior.at(k).x, expected_x, EXP_NEAR_TOLERANCE)
        << "x mismatch at interior node " << k;
    EXPECT_NEAR(interior.at(k).y, expected_y, EXP_NEAR_TOLERANCE)
        << "y mismatch at interior node " << k;
    EXPECT_NEAR(interior.at(k).z, expected_z, EXP_NEAR_TOLERANCE)
        << "z mismatch at interior node " << k;
  }
}

TEST_F(SegmentTest, MidpointQuadratureAbscissaMapsToSegmentMidpoint) {
  // The 3-point GLL rule is exactly {-1, 0, 1}, so its single interior point
  // (xi = 0) should map to the geometric midpoint of the segment -- this
  // follows directly from any correct reference-to-physical mapping,
  // independent of the exact formula used.
  const Node left = MakeNode(0, 0.0, 0.0, 0.0);
  const Node right = MakeNode(1, 10.0, 20.0, -6.0);
  std::vector<Node> nodes = {left, right};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 1, 0, mesh);
  GaussLobattoLegendre gll(3);

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_EQ(interior.size(), 1u);
  EXPECT_NEAR(interior.at(0).x, 5.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(interior.at(0).y, 10.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(interior.at(0).z, -3.0, EXP_NEAR_TOLERANCE);
}

TEST_F(SegmentTest, InteriorNodesAreMonotonicAlongTheSegmentDirection) {
  // Regardless of the exact mapping formula, successive interior nodes
  // (which come from successive, increasing GLL abscissas) should move
  // monotonically from the left endpoint toward the right one.
  const Node left = MakeNode(0, 0.0, 0.0, 0.0);
  const Node right = MakeNode(1, 10.0, 0.0, 0.0);
  std::vector<Node> nodes = {left, right};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 1, 0, mesh);
  GaussLobattoLegendre gll(6);

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_GE(interior.size(), 2u);
  for (size_t k = 0; k + 1 < interior.size(); ++k) {
    EXPECT_LT(interior.at(k).x, interior.at(k + 1).x)
        << "Interior nodes are not monotonically increasing along x.";
  }
}

TEST_F(SegmentTest, DirectionFollowsBCNodeIdOrder) {
  // Swapping which ID is boundary_node_ids_[0] vs [1] should reverse which
  // endpoint is "left" for the purposes of the mapping.
  const Node node_a = MakeNode(0, 0.0, 0.0, 0.0);
  const Node node_b = MakeNode(1, 10.0, 0.0, 0.0);
  std::vector<Node> nodes = {node_a, node_b};
  mesh.AddNodes(nodes);
  Segment segment({1, 0}, 1, 0, mesh);  // node_b is "left", node_a is "right"
  GaussLobattoLegendre gll(3);

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_EQ(interior.size(), 1u);
  // The midpoint is the same regardless of direction, so this alone can't
  // distinguish a correct implementation from a reversed one -- see
  // InteriorNodesAreMonotonicAlongTheSegmentDirection-style checks with a
  // non-midpoint quadrature for that.
  EXPECT_NEAR(interior.at(0).x, 5.0, EXP_NEAR_TOLERANCE);
}

// ---------------------------------------------------------------------------
// BC inheritance
// ---------------------------------------------------------------------------

TEST_F(SegmentTest, InteriorBCMatchesSharedEndpointBC) {
  const unsigned int kSharedBCID = 7u;
  Node left = MakeNode(0, 0.0, 0.0, 0.0, kSharedBCID);
  Node right = MakeNode(1, 10.0, 0.0, 0.0, kSharedBCID);
  std::vector<Node> nodes = {left, right};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 1, 0, mesh);
  GaussLobattoLegendre gll(5);  // 3 interior points

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_FALSE(interior.empty());
  for (const auto& node : interior) {
    EXPECT_EQ(node.bc_id, kSharedBCID);
  }
}

TEST_F(SegmentTest, InteriorBCIsNoneWhenEndpointBoundariesDiffer) {
  Node left = MakeNode(0, 0.0, 0.0, 0.0, 0u);
  Node right = MakeNode(1, 10.0, 0.0, 0.0, 7u);
  std::vector<Node> nodes = {left, right};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 1, 0, mesh);
  GaussLobattoLegendre gll(5);  // 3 interior points

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_FALSE(interior.empty());
  for (const auto& node : interior) {
    EXPECT_EQ(node.bc_id, 0u);
  }
}

TEST_F(SegmentTest, InteriorBCIsNoneWhenBothEndpointsAreNone) {
  // Degenerate but common case: neither endpoint is on a boundary, so
  // interior nodes shouldn't be either.
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0, 0u),
                             MakeNode(1, 10.0, 0.0, 0.0, 0u)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 1, 0, mesh);
  GaussLobattoLegendre gll(4);  // 2 interior points

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_FALSE(interior.empty());
  for (const auto& node : interior) {
    EXPECT_EQ(node.bc_id, 0u);
  }
}

// ---------------------------------------------------------------------------
// LocalStiffnessMatrix
//
// Per Segment::LocalStiffnessMatrix's documented formula,
//   K_ij = (Omega_x^2 / Sigma_t) * (2 / h_e) * sum_k w_k * dl_i/dxi(xi_k) *
//          dl_j/dxi(xi_k)
// ---------------------------------------------------------------------------

class SegmentStiffnessMatrixTest : public testing::Test {
 protected:
  Mesh mesh;
};

TEST_F(SegmentStiffnessMatrixTest, TwoPointLinearElementMatchesHandDerivedValues) {
  // N=2: dl_0/dxi = -0.5 and dl_1/dxi = 0.5 everywhere (linear basis), and
  // the 2-point GLL rule has unit weights. This is the textbook 1D linear
  // finite-element stiffness matrix:
  //   K = (Omega_x^2 / Sigma_t) * (1 / h_e) * [[1, -1], [-1, 1]]
  const double length = 4.0;
  const double sigma_t = 2.0;
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, length, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 0, 0, mesh);
  GaussLobattoLegendre gll(2);
  MaterialBank material_bank = MakeSingleMaterialBank(sigma_t);
  Ordinate ordinate = MakeOrdinateWithXCosine(1.0);

  auto stiffness = segment.LocalStiffnessMatrix(gll, material_bank, ordinate);

  const double expected_magnitude = (1.0 / sigma_t) * (1.0 / length);
  EXPECT_NEAR(stiffness(0, 0), expected_magnitude, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(stiffness(1, 1), expected_magnitude, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(stiffness(0, 1), -expected_magnitude, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(stiffness(1, 0), -expected_magnitude, EXP_NEAR_TOLERANCE);
}

TEST_F(SegmentStiffnessMatrixTest,
       ThreePointQuadraticElementMatchesHandDerivedValues) {
  // N=3 (quadratic spectral element): with abscissas {-1, 0, 1}, weights
  // {1/3, 4/3, 1/3}, and the well-known closed-form Lagrange derivatives at
  // those nodes, the correct GLL(3) stiffness matrix (before the
  // Omega_x^2/Sigma_t * 2/h_e prefactor) is the classic
  //   [[ 7/6, -4/3,  1/6],
  //    [-4/3,  8/3, -4/3],
  //    [ 1/6, -4/3,  7/6]]
  const double length = 4.0;
  const double sigma_t = 2.0;
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, length, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 0, 0, mesh);
  GaussLobattoLegendre gll(3);
  MaterialBank material_bank = MakeSingleMaterialBank(sigma_t);
  Ordinate ordinate = MakeOrdinateWithXCosine(1.0);

  auto stiffness = segment.LocalStiffnessMatrix(gll, material_bank, ordinate);

  const double coeff = (1.0 / sigma_t) * (2.0 / length);
  EXPECT_NEAR(stiffness(0, 0), coeff * 7.0 / 6.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(stiffness(0, 1), coeff * -4.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(stiffness(0, 2), coeff * 1.0 / 6.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(stiffness(1, 0), coeff * -4.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(stiffness(1, 1), coeff * 8.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(stiffness(1, 2), coeff * -4.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(stiffness(2, 0), coeff * 1.0 / 6.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(stiffness(2, 1), coeff * -4.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(stiffness(2, 2), coeff * 7.0 / 6.0, EXP_NEAR_TOLERANCE);
}

class SegmentStiffnessMatrixOrderTest
    : public ::testing::TestWithParam<size_t> {
 protected:
  Mesh mesh;
};

TEST_P(SegmentStiffnessMatrixOrderTest, IsSymmetric) {
  const size_t n_points = GetParam();
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 3.0, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 0, 0, mesh);
  GaussLobattoLegendre gll(n_points);
  MaterialBank material_bank = MakeSingleMaterialBank(1.5);
  Ordinate ordinate = MakeOrdinateWithXCosine(0.8);

  auto stiffness = segment.LocalStiffnessMatrix(gll, material_bank, ordinate);

  for (size_t i = 0; i < n_points; ++i)
    for (size_t j = 0; j < n_points; ++j)
      EXPECT_DOUBLE_EQ(stiffness(i, j), stiffness(j, i))
          << "n=" << n_points << ", (i,j)=(" << i << "," << j << ")";
}

TEST_P(SegmentStiffnessMatrixOrderTest, RowsSumToZero) {
  // Physically, K applied to the constant vector must be zero: the
  // derivative of a constant is zero, so the weak-form stiffness operator
  // annihilates constants regardless of element order.
  const size_t n_points = GetParam();
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 3.0, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 0, 0, mesh);
  GaussLobattoLegendre gll(n_points);
  MaterialBank material_bank = MakeSingleMaterialBank(1.5);
  Ordinate ordinate = MakeOrdinateWithXCosine(0.8);

  auto stiffness = segment.LocalStiffnessMatrix(gll, material_bank, ordinate);

  for (size_t i = 0; i < n_points; ++i) {
    double row_sum = 0.0;
    for (size_t j = 0; j < n_points; ++j) row_sum += stiffness(i, j);
    EXPECT_NEAR(row_sum, 0.0, 1e-10) << "n=" << n_points << ", row " << i;
  }
}

INSTANTIATE_TEST_SUITE_P(VariousOrders, SegmentStiffnessMatrixOrderTest,
                         ::testing::Values(2, 3, 4, 5, 6, 8));

TEST_F(SegmentStiffnessMatrixTest, ScalesInverselyWithLength) {
  MaterialBank material_bank = MakeSingleMaterialBank(1.0);
  Ordinate ordinate = MakeOrdinateWithXCosine(1.0);
  GaussLobattoLegendre gll(4);

  Mesh short_mesh;
  short_mesh.AddNodes({MakeNode(0, 0.0, 0.0, 0.0), MakeNode(1, 2.0, 0.0, 0.0)});
  Segment short_segment({0, 1}, 0, 0, short_mesh);

  Mesh long_mesh;
  long_mesh.AddNodes({MakeNode(0, 0.0, 0.0, 0.0), MakeNode(1, 4.0, 0.0, 0.0)});
  Segment long_segment({0, 1}, 0, 0, long_mesh);

  auto k_short = short_segment.LocalStiffnessMatrix(gll, material_bank, ordinate);
  auto k_long = long_segment.LocalStiffnessMatrix(gll, material_bank, ordinate);

  for (size_t i = 0; i < gll.n_points(); ++i)
    for (size_t j = 0; j < gll.n_points(); ++j)
      EXPECT_NEAR(k_long(i, j), k_short(i, j) / 2.0, 1e-12);
}

TEST_F(SegmentStiffnessMatrixTest, ScalesInverselyWithTotalXS) {
  Ordinate ordinate = MakeOrdinateWithXCosine(1.0);
  GaussLobattoLegendre gll(4);
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 5.0, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 0, 0, mesh);

  MaterialBank bank_1x = MakeSingleMaterialBank(1.0);
  MaterialBank bank_2x = MakeSingleMaterialBank(2.0);

  auto k_1x = segment.LocalStiffnessMatrix(gll, bank_1x, ordinate);
  auto k_2x = segment.LocalStiffnessMatrix(gll, bank_2x, ordinate);

  for (size_t i = 0; i < gll.n_points(); ++i)
    for (size_t j = 0; j < gll.n_points(); ++j)
      EXPECT_NEAR(k_2x(i, j), k_1x(i, j) / 2.0, 1e-12);
}

TEST_F(SegmentStiffnessMatrixTest, ScalesWithOmegaXSquared) {
  MaterialBank material_bank = MakeSingleMaterialBank(1.0);
  GaussLobattoLegendre gll(4);
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 5.0, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 0, 0, mesh);

  Ordinate full = MakeOrdinateWithXCosine(1.0);
  Ordinate half = MakeOrdinateWithXCosine(0.5);

  auto k_full = segment.LocalStiffnessMatrix(gll, material_bank, full);
  auto k_half = segment.LocalStiffnessMatrix(gll, material_bank, half);

  for (size_t i = 0; i < gll.n_points(); ++i)
    for (size_t j = 0; j < gll.n_points(); ++j)
      EXPECT_NEAR(k_half(i, j), k_full(i, j) * 0.25, 1e-12);
}

// ---------------------------------------------------------------------------
// LocalMassMatrix
//
// Per Segment::LocalMassMatrix's documented formula, the mass matrix is
// diagonal with M_ii = Sigma_t * (h_e / 2) * w_i.
// ---------------------------------------------------------------------------

class SegmentMassMatrixTest : public testing::Test {
 protected:
  Mesh mesh;
};

TEST_F(SegmentMassMatrixTest, IsDiagonal) {
  const double length = 3.0;
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, length, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 0, 0, mesh);
  GaussLobattoLegendre gll(5);
  MaterialBank material_bank = MakeSingleMaterialBank(1.0);

  auto mass = segment.LocalMassMatrix(gll, material_bank);

  for (size_t i = 0; i < gll.n_points(); ++i)
    for (size_t j = 0; j < gll.n_points(); ++j)
      if (i != j) EXPECT_DOUBLE_EQ(mass(i, j), 0.0) << "(i,j)=(" << i << "," << j << ")";
}

TEST_F(SegmentMassMatrixTest, TwoPointMatchesAnalyticalValue) {
  // N=2: both GLL weights are exactly 1, so M = Sigma_t * h_e / 2 * I.
  const double length = 5.0;
  const double sigma_t = 3.0;
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, length, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 0, 0, mesh);
  GaussLobattoLegendre gll(2);
  MaterialBank material_bank = MakeSingleMaterialBank(sigma_t);

  auto mass = segment.LocalMassMatrix(gll, material_bank);

  const double expected = sigma_t * length / 2.0;
  EXPECT_DOUBLE_EQ(mass(0, 0), expected);
  EXPECT_DOUBLE_EQ(mass(1, 1), expected);
}

TEST_F(SegmentMassMatrixTest, ThreePointMatchesHandDerivedWeights) {
  // N=3: GLL weights are exactly {1/3, 4/3, 1/3}.
  const double length = 6.0;
  const double sigma_t = 2.0;
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, length, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 0, 0, mesh);
  GaussLobattoLegendre gll(3);
  MaterialBank material_bank = MakeSingleMaterialBank(sigma_t);

  auto mass = segment.LocalMassMatrix(gll, material_bank);

  const double half_sigma_h = sigma_t * length / 2.0;
  EXPECT_NEAR(mass(0, 0), half_sigma_h * (1.0 / 3.0), EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(mass(1, 1), half_sigma_h * (4.0 / 3.0), EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(mass(2, 2), half_sigma_h * (1.0 / 3.0), EXP_NEAR_TOLERANCE);
}

class SegmentMassMatrixOrderTest : public ::testing::TestWithParam<size_t> {
 protected:
  Mesh mesh;
};

TEST_P(SegmentMassMatrixOrderTest, TraceEqualsTotalXsTimesLength) {
  // trace(M) = Sigma_t * h_e / 2 * sum_i(w_i) = Sigma_t * h_e / 2 * 2 =
  // Sigma_t * h_e, since GLL weights always sum to exactly 2 regardless of
  // order.
  const size_t n_points = GetParam();
  const double length = 7.0;
  const double sigma_t = 4.0;
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, length, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 0, 0, mesh);
  GaussLobattoLegendre gll(n_points);
  MaterialBank material_bank = MakeSingleMaterialBank(sigma_t);

  auto mass = segment.LocalMassMatrix(gll, material_bank);

  double trace = 0.0;
  for (size_t i = 0; i < n_points; ++i) trace += mass(i, i);
  EXPECT_NEAR(trace, sigma_t * length, 1e-12) << "n=" << n_points;
}

INSTANTIATE_TEST_SUITE_P(VariousOrders, SegmentMassMatrixOrderTest,
                         ::testing::Values(2, 3, 4, 5, 6, 8));

TEST_F(SegmentMassMatrixTest, ScalesWithLength) {
  MaterialBank material_bank = MakeSingleMaterialBank(1.0);
  GaussLobattoLegendre gll(4);

  Mesh short_mesh;
  short_mesh.AddNodes({MakeNode(0, 0.0, 0.0, 0.0), MakeNode(1, 2.0, 0.0, 0.0)});
  Segment short_segment({0, 1}, 0, 0, short_mesh);

  Mesh long_mesh;
  long_mesh.AddNodes({MakeNode(0, 0.0, 0.0, 0.0), MakeNode(1, 4.0, 0.0, 0.0)});
  Segment long_segment({0, 1}, 0, 0, long_mesh);

  auto m_short = short_segment.LocalMassMatrix(gll, material_bank);
  auto m_long = long_segment.LocalMassMatrix(gll, material_bank);

  for (size_t i = 0; i < gll.n_points(); ++i)
    EXPECT_NEAR(m_long(i, i), m_short(i, i) * 2.0, 1e-12);
}

TEST_F(SegmentMassMatrixTest, ScalesWithTotalXS) {
  GaussLobattoLegendre gll(4);
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 5.0, 0.0, 0.0)};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 0, 0, mesh);

  MaterialBank bank_1x = MakeSingleMaterialBank(1.0);
  MaterialBank bank_3x = MakeSingleMaterialBank(3.0);

  auto m_1x = segment.LocalMassMatrix(gll, bank_1x);
  auto m_3x = segment.LocalMassMatrix(gll, bank_3x);

  for (size_t i = 0; i < gll.n_points(); ++i)
    EXPECT_NEAR(m_3x(i, i), m_1x(i, i) * 3.0, 1e-12);
}

// ---------------------------------------------------------------------------
// LocalForcingVector
//
// Per Segment::LocalForcingVector's documented formula,
//   f_i = h_e/2 * sum_k w_k * Q(xi_k) * l_i(xi_k)
//         + sum_k w_k * Q(xi_k) * dl_i/dxi(xi_k)
// which, by the cardinality of the Lagrange basis (l_i(xi_k) = delta_ik),
// simplifies the first term to h_e/2 * w_i * Q(xi_i).
//
// These tests call CreateInteriorNodes through the real public API to
// populate node_ids_ (see SegmentTest.CreateInteriorNodesPopulatesNodeIds),
// rather than reaching into Segment's internals.
// ---------------------------------------------------------------------------

class SegmentForcingVectorTest : public testing::Test {
 protected:
  Mesh mesh;
};

TEST_F(SegmentForcingVectorTest, TwoPointDistinctBoundarySourcesMatchHandDerivedValues) {
  // N=2: dl_0/dxi = -0.5 and dl_1/dxi = 0.5 everywhere, weights are both 1.
  //   f_0 = h/2*Q_0 + (-0.5)*(Q_0 + Q_1)
  //   f_1 = h/2*Q_1 + ( 0.5)*(Q_0 + Q_1)
  // With h=4, Q_0=2, Q_1=6: f_0 = 4 - 4 = 0, f_1 = 12 + 4 = 16.
  const double length = 4.0;
  const double q0 = 2.0;
  const double q1 = 6.0;
  Node left = MakeNode(0, 0.0, 0.0, 0.0);
  left.source_fluxes = {q0};
  Node right = MakeNode(1, length, 0.0, 0.0);
  right.source_fluxes = {q1};
  mesh.AddNodes({left, right});
  Segment segment({0, 1}, 0, 0, mesh);
  GaussLobattoLegendre gll(2);
  segment.CreateInteriorNodes(mesh.nodes(), gll);  // no interior points; N=2

  auto forcing = segment.LocalForcingVector(gll, mesh, 0);

  EXPECT_NEAR(forcing(0), 0.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(forcing(1), 16.0, EXP_NEAR_TOLERANCE);
}

TEST_F(SegmentForcingVectorTest, TwoPointZeroSourceProducesZeroVector) {
  Node left = MakeNode(0, 0.0, 0.0, 0.0);
  left.source_fluxes = {0.0};
  Node right = MakeNode(1, 4.0, 0.0, 0.0);
  right.source_fluxes = {0.0};
  mesh.AddNodes({left, right});
  Segment segment({0, 1}, 0, 0, mesh);
  GaussLobattoLegendre gll(2);
  segment.CreateInteriorNodes(mesh.nodes(), gll);  // no interior points; N=2

  auto forcing = segment.LocalForcingVector(gll, mesh, 0);

  EXPECT_DOUBLE_EQ(forcing(0), 0.0);
  EXPECT_DOUBLE_EQ(forcing(1), 0.0);
}

TEST_F(SegmentForcingVectorTest,
       ThreePointSingleBoundarySourceMatchesHandDerivedValues) {
  // N=3, source Q = [0, 0, 1] (only the right boundary node emits). Only the
  // k=2 quadrature term survives:
  //   f_i = [i==2] * h/2 * w_2 + w_2 * dl_i/dxi(xi_2)
  // With h=6, w_2=1/3, and dl_i/dxi(1) = {0.5, -2.0, 1.5} for i={0,1,2}:
  //   f_0 = 1/3 * 0.5        = 1/6
  //   f_1 = 1/3 * -2.0       = -2/3
  //   f_2 = 6/2*1/3 + 1/3*1.5 = 1 + 0.5 = 1.5
  const double length = 6.0;
  Node left = MakeNode(0, 0.0, 0.0, 0.0);
  left.source_fluxes = {0.0};
  Node right = MakeNode(1, length, 0.0, 0.0);
  right.source_fluxes = {1.0};
  mesh.AddNodes({left, right});
  Segment segment({0, 1}, 0, 0, mesh);
  GaussLobattoLegendre gll(3);
  auto interior = segment.CreateInteriorNodes(mesh.nodes(), gll);
  ASSERT_EQ(interior.size(), 1u);
  interior.at(0).source_fluxes = {0.0};
  mesh.AddNodes(interior);

  auto forcing = segment.LocalForcingVector(gll, mesh, 0);

  EXPECT_NEAR(forcing(0), 1.0 / 6.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(forcing(1), -2.0 / 3.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(forcing(2), 1.5, EXP_NEAR_TOLERANCE);
}

class SegmentForcingVectorConstantSourceTest
    : public ::testing::TestWithParam<size_t> {
 protected:
  Mesh mesh;
};

TEST_P(SegmentForcingVectorConstantSourceTest,
       MatchesClosedFormFormulaAcrossOrders) {
  // For a spatially-constant source Q0, the derivative term integrates
  // exactly (GLL is exact for the low-degree dl_i/dxi):
  //   sum_k w_k * dl_i/dxi(xi_k) = l_i(1) - l_i(-1) = [i==N-1] - [i==0]
  // so f_i = h/2 * w_i * Q0 + Q0 * ([i==N-1] - [i==0]).
  const size_t n_points = GetParam();
  const double length = 5.0;
  const double q0 = 3.0;

  GaussLobattoLegendre gll(n_points);
  Node left = MakeNode(0, 0.0, 0.0, 0.0);
  left.source_fluxes = {q0};
  Node right = MakeNode(1, length, 0.0, 0.0);
  right.source_fluxes = {q0};
  mesh.AddNodes({left, right});
  Segment segment({0, 1}, 0, 0, mesh);
  auto interior = segment.CreateInteriorNodes(mesh.nodes(), gll);
  for (auto& node : interior) node.source_fluxes = {q0};
  mesh.AddNodes(interior);

  auto forcing = segment.LocalForcingVector(gll, mesh, 0);

  for (size_t i = 0; i < n_points; ++i) {
    double boundary_term = 0.0;
    if (i == n_points - 1) boundary_term = q0;
    if (i == 0) boundary_term -= q0;
    const double expected = length / 2.0 * gll.GetWeight(i) * q0 + boundary_term;
    EXPECT_NEAR(forcing(i), expected, 1e-10) << "n=" << n_points << ", i=" << i;
  }
}

INSTANTIATE_TEST_SUITE_P(VariousOrders,
                         SegmentForcingVectorConstantSourceTest,
                         ::testing::Values(2, 3, 4, 5, 6));

TEST_F(SegmentForcingVectorTest, IsLinearInSourceValues) {
  // f(Q_a + Q_b) == f(Q_a) + f(Q_b), since every term in the documented
  // formula is linear in the nodal source values. sources.at(k) is the
  // source at local GLL position k (left, interior..., right).
  const double length = 5.0;
  GaussLobattoLegendre gll(4);

  auto build_and_evaluate = [&](std::vector<double> sources) {
    Mesh local_mesh;
    Node left = MakeNode(0, 0.0, 0.0, 0.0);
    left.source_fluxes = {sources.at(0)};
    Node right = MakeNode(1, length, 0.0, 0.0);
    right.source_fluxes = {sources.back()};
    local_mesh.AddNodes({left, right});
    Segment segment({0, 1}, 0, 0, local_mesh);
    auto interior = segment.CreateInteriorNodes(local_mesh.nodes(), gll);
    for (size_t k = 0; k < interior.size(); ++k)
      interior.at(k).source_fluxes = {sources.at(k + 1)};
    local_mesh.AddNodes(interior);
    return segment.LocalForcingVector(gll, local_mesh, 0);
  };

  std::vector<double> q_a = {1.0, -2.0, 3.0, 0.5};
  std::vector<double> q_b = {4.0, 1.5, -1.0, 2.0};
  std::vector<double> q_sum = {5.0, -0.5, 2.0, 2.5};

  auto f_a = build_and_evaluate(q_a);
  auto f_b = build_and_evaluate(q_b);
  auto f_sum = build_and_evaluate(q_sum);

  for (size_t i = 0; i < gll.n_points(); ++i)
    EXPECT_NEAR(f_sum(i), f_a(i) + f_b(i), 1e-10) << "i=" << i;
}

TEST_F(SegmentForcingVectorTest,
       LocalForcingVectorThrowsWhenCalledBeforeCreateInteriorNodes) {
  // node_ids_ is populated as a side effect of CreateInteriorNodes (see
  // SegmentTest.CreateInteriorNodesPopulatesNodeIds), so a Segment built but
  // not yet prepared has an empty node_ids_ and LocalForcingVector correctly
  // fails rather than silently reading garbage.
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 4.0, 0.0, 0.0)};
  nodes.at(0).source_fluxes = {1.0};
  nodes.at(1).source_fluxes = {1.0};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 0, 0, mesh);
  GaussLobattoLegendre gll(2);

  EXPECT_TRUE(segment.node_ids().empty());
  EXPECT_THROW(segment.LocalForcingVector(gll, mesh, 0), std::out_of_range);
}

TEST_F(SegmentForcingVectorTest,
       SucceedsAfterCreateInteriorNodesPopulatesNodeIdsThroughPublicAPI) {
  // Regression test: once CreateInteriorNodes has populated node_ids_
  // through the ordinary public API (no test-only backdoor), a Segment can
  // successfully compute its local forcing vector, stiffness matrix, and
  // mass matrix.
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 4.0, 0.0, 0.0)};
  nodes.at(0).source_fluxes = {1.0};
  nodes.at(1).source_fluxes = {1.0};
  mesh.AddNodes(nodes);
  Segment segment({0, 1}, 0, 0, mesh);
  GaussLobattoLegendre gll(3);

  auto interior = segment.CreateInteriorNodes(mesh.nodes(), gll);
  for (auto& node : interior) node.source_fluxes = {1.0};
  mesh.AddNodes(interior);

  EXPECT_FALSE(segment.node_ids().empty());
  EXPECT_NO_THROW(segment.LocalForcingVector(gll, mesh, 0));

  MaterialBank material_bank = MakeSingleMaterialBank(1.0);
  Ordinate ordinate = MakeOrdinateWithXCosine(1.0);
  EXPECT_NO_THROW(segment.LocalStiffnessMatrix(gll, material_bank, ordinate));
  EXPECT_NO_THROW(segment.LocalMassMatrix(gll, material_bank));
}

}  // namespace hummingbird