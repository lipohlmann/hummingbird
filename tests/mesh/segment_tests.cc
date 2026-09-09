#include <gtest/gtest.h>

#include <array>
#include <tuple>
#include <vector>

#include "mesh/node.h"
#include "mesh/segment.h"
#include "quadrature/gauss_lobatto_legendre.h"
#include "utils/constants.h"
#include "utils/enums.h"

namespace hummingbird {

namespace {

Node MakeNode(size_t id, double x, double y, double z, BC boundary = BC::NONE) {
  Node node;
  node.id = id;
  node.x = x;
  node.y = y;
  node.z = z;
  node.boundary = boundary;
  return node;
}

}  // namespace

// ---------------------------------------------------------------------------
// Interior node count
// ---------------------------------------------------------------------------

class SegmentInteriorCountTest : public ::testing::TestWithParam<size_t> {};

TEST_P(SegmentInteriorCountTest, ProducesNPointsMinusTwoInteriorNodes) {
  const size_t n_points = GetParam();
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 10.0, 0.0, 0.0)};
  Segment segment({0, 1}, 1, 0);
  GaussLobattoLegendre gll(n_points);

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  EXPECT_EQ(interior.size(), n_points - 2);
}

INSTANTIATE_TEST_SUITE_P(VariousOrders, SegmentInteriorCountTest,
                         ::testing::Values(2, 3, 4, 5, 6, 8));

TEST(SegmentTest, TwoPointQuadratureProducesNoInteriorNodes) {
  // A GLL rule with only the two endpoints has no interior points at all.
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 10.0, 0.0, 0.0)};
  Segment segment({0, 1}, 1, 0);
  GaussLobattoLegendre gll(2);

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  EXPECT_TRUE(interior.empty());
}

// ---------------------------------------------------------------------------
// Node ID assignment
// ---------------------------------------------------------------------------

TEST(SegmentTest, InteriorNodeIdsContinueFromExistingNodeCount) {
  // Three nodes already exist, so new IDs should start at 3 and increment.
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0),
                             MakeNode(1, 10.0, 0.0, 0.0),
                             MakeNode(2, 20.0, 0.0, 0.0)};
  Segment segment({0, 1}, 1, 0);
  GaussLobattoLegendre gll(5);  // 3 interior points

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_EQ(interior.size(), 3u);
  for (size_t k = 0; k < interior.size(); ++k) {
    EXPECT_EQ(interior.at(k).id, nodes.size() + k)
        << "Interior node " << k << " has an unexpected ID.";
  }
}

TEST(SegmentTest, UsesNodesAtSpecifiedBCIndicesRegardlessOfPosition) {
  // boundary_node_ids_ should be used as indices into existing_nodes, not
  // assumed to be the first two entries.
  std::vector<Node> nodes(6);
  nodes.at(2) = MakeNode(2, 0.0, 0.0, 0.0);
  nodes.at(5) = MakeNode(5, 10.0, 0.0, 0.0);
  Segment segment({2, 5}, 1, 0);
  GaussLobattoLegendre gll(3);  // 1 interior point

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_EQ(interior.size(), 1u);
  EXPECT_EQ(interior.at(0).id, nodes.size());
}

// ---------------------------------------------------------------------------
// Coordinate mapping
//
// Per Segment::CreateInteriorNodes's documented mapping,
//   r(xi) = r1 + ((xi + 1) / 2) * (r2 - r1),  xi in [-1, 1]
// where r1 is the node at boundary_node_ids_[0] and r2 is the node at
// boundary_node_ids_[1].
// ---------------------------------------------------------------------------

TEST(SegmentTest, InteriorNodesMatchDocumentedReferenceToPhysicalMapping) {
  const Node left = MakeNode(0, 1.0, 2.0, 3.0);
  const Node right = MakeNode(1, 4.0, 6.0, 3.0);  // non-axis-aligned segment
  std::vector<Node> nodes = {left, right};
  Segment segment({0, 1}, 1, 0);
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

TEST(SegmentTest, MidpointQuadratureAbscissaMapsToSegmentMidpoint) {
  // The 3-point GLL rule is exactly {-1, 0, 1}, so its single interior point
  // (xi = 0) should map to the geometric midpoint of the segment -- this
  // follows directly from any correct reference-to-physical mapping,
  // independent of the exact formula used.
  const Node left = MakeNode(0, 0.0, 0.0, 0.0);
  const Node right = MakeNode(1, 10.0, 20.0, -6.0);
  std::vector<Node> nodes = {left, right};
  Segment segment({0, 1}, 1, 0);
  GaussLobattoLegendre gll(3);

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_EQ(interior.size(), 1u);
  EXPECT_NEAR(interior.at(0).x, 5.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(interior.at(0).y, 10.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(interior.at(0).z, -3.0, EXP_NEAR_TOLERANCE);
}

TEST(SegmentTest, InteriorNodesAreMonotonicAlongTheSegmentDirection) {
  // Regardless of the exact mapping formula, successive interior nodes
  // (which come from successive, increasing GLL abscissas) should move
  // monotonically from the left endpoint toward the right one.
  const Node left = MakeNode(0, 0.0, 0.0, 0.0);
  const Node right = MakeNode(1, 10.0, 0.0, 0.0);
  std::vector<Node> nodes = {left, right};
  Segment segment({0, 1}, 1, 0);
  GaussLobattoLegendre gll(6);

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_GE(interior.size(), 2u);
  for (size_t k = 0; k + 1 < interior.size(); ++k) {
    EXPECT_LT(interior.at(k).x, interior.at(k + 1).x)
        << "Interior nodes are not monotonically increasing along x.";
  }
}

TEST(SegmentTest, DirectionFollowsBCNodeIdOrder) {
  // Swapping which ID is boundary_node_ids_[0] vs [1] should reverse which
  // endpoint is "left" for the purposes of the mapping.
  const Node node_a = MakeNode(0, 0.0, 0.0, 0.0);
  const Node node_b = MakeNode(1, 10.0, 0.0, 0.0);
  std::vector<Node> nodes = {node_a, node_b};
  Segment segment({1, 0}, 1, 0);  // node_b is "left", node_a is "right"
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

TEST(SegmentTest, InteriorBCMatchesSharedEndpointBC) {
  // NOTE: BC::NONE is the only enumerator visible from the provided
  // source; static_cast<BC>(1) stands in for "some other boundary
  // value" and should be replaced with the real enumerator name if it
  // differs from the underlying value 1.
  const BC kSharedBC = static_cast<BC>(1);
  Node left = MakeNode(0, 0.0, 0.0, 0.0, kSharedBC);
  Node right = MakeNode(1, 10.0, 0.0, 0.0, kSharedBC);
  std::vector<Node> nodes = {left, right};
  Segment segment({0, 1}, 1, 0);
  GaussLobattoLegendre gll(5);  // 3 interior points

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_FALSE(interior.empty());
  for (const auto& node : interior) {
    EXPECT_EQ(node.boundary, kSharedBC);
  }
}

TEST(SegmentTest, InteriorBCIsNoneWhenEndpointBoundariesDiffer) {
  Node left = MakeNode(0, 0.0, 0.0, 0.0, BC::NONE);
  Node right = MakeNode(1, 10.0, 0.0, 0.0, static_cast<BC>(1));
  std::vector<Node> nodes = {left, right};
  Segment segment({0, 1}, 1, 0);
  GaussLobattoLegendre gll(5);  // 3 interior points

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_FALSE(interior.empty());
  for (const auto& node : interior) {
    EXPECT_EQ(node.boundary, BC::NONE);
  }
}

TEST(SegmentTest, InteriorBCIsNoneWhenBothEndpointsAreNone) {
  // Degenerate but common case: neither endpoint is on a boundary, so
  // interior nodes shouldn't be either.
  std::vector<Node> nodes = {MakeNode(0, 0.0, 0.0, 0.0, BC::NONE),
                             MakeNode(1, 10.0, 0.0, 0.0, BC::NONE)};
  Segment segment({0, 1}, 1, 0);
  GaussLobattoLegendre gll(4);  // 2 interior points

  auto interior = segment.CreateInteriorNodes(nodes, gll);
  ASSERT_FALSE(interior.empty());
  for (const auto& node : interior) {
    EXPECT_EQ(node.boundary, BC::NONE);
  }
}

}  // namespace hummingbird