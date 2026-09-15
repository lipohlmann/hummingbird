#include "mesh/mesh.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

#include "mesh/segment.h"
#include "quadrature/gauss_lobatto_legendre.h"
#include "utils/constants.h"

namespace hummingbird {

namespace {

// Mirrors cases/1D/1D_MMS_1_gmsh.msh: a unit-length 1D domain split into 4
// segments, with "bc_west"/"bc_east" boundary points and a single
// "material:mms_material" / "source:mms_source" curve region.
constexpr char kOneDGmsh[] = R"(
$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
4
0 2 "bc_west:vacuum"
0 3 "bc_east:vacuum"
1 4 "material:mms_material"
1 5 "source:mms_source"
$EndPhysicalNames
$Entities
2 1 0 0
1 0 0 0 1 2
2 1 0 0 1 3
1 -9.999999994736442e-08 -1e-07 -1e-07 1.0000001 1e-07 1e-07 2 4 5 2 2 -1
$EndEntities
$Nodes
3 5 1 5
0 1 0 1
1
0 0 0
0 2 0 1
2
1 0 0
1 1 0 3
3
4
5
0.75 0 0
0.5 0 0
0.25 0 0
$EndNodes
$Elements
3 6 1 6
0 1 15 1
1 1
0 2 15 1
2 2
1 1 1 4
3 2 3
4 3 4
5 4 5
6 5 1
$EndElements
)";

std::string WriteTempMesh(const std::string& name, const std::string& contents) {
  std::string path = ::testing::TempDir() + name;
  std::ofstream file(path);
  file << contents;
  return path;
}

Node MakeNode(size_t id, double x, double y, double z) {
  Node node;
  node.id = id;
  node.x = x;
  node.y = y;
  node.z = z;
  return node;
}

}  // namespace

TEST(MeshGMSHTest, ReadsExpectedNumberOfNodesAndElements) {
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_basic.msh", kOneDGmsh));
  EXPECT_EQ(mesh.nodes().size(), 5u);
  EXPECT_EQ(mesh.n_elements(), 4u);
}

TEST(MeshGMSHTest, NodeCoordinatesMatchMeshFile) {
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_coords.msh", kOneDGmsh));
  const auto& nodes = mesh.nodes();
  ASSERT_EQ(nodes.size(), 5u);
  EXPECT_NEAR(nodes.at(0).x, 0.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(nodes.at(1).x, 1.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(nodes.at(2).x, 0.75, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(nodes.at(3).x, 0.5, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(nodes.at(4).x, 0.25, EXP_NEAR_TOLERANCE);
}

TEST(MeshGMSHTest, ElementsHaveMaterialIDFromMaterialPhysicalGroup) {
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_material.msh", kOneDGmsh));
  ASSERT_EQ(mesh.n_elements(), 4u);
  for (size_t i = 0; i < mesh.n_elements(); i++)
    EXPECT_EQ(mesh.GetElement(i).material_id(), 4)
        << "Element " << i << " has an unexpected material ID.";
}

TEST(MeshGMSHTest, ElementsHaveSourceIDFromSourcePhysicalGroup) {
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_source.msh", kOneDGmsh));
  ASSERT_EQ(mesh.n_elements(), 4u);
  for (size_t i = 0; i < mesh.n_elements(); i++)
    EXPECT_EQ(mesh.GetElement(i).source_id(), 5)
        << "Element " << i << " has an unexpected source ID.";
}

TEST(MeshGMSHTest, SegmentConnectivityMatchesFileOrder) {
  // Segment::CreateInteriorNodes maps each element's boundary nodes to the
  // segment midpoint when using a 3-point GLL rule, so the resulting
  // interior node coordinates confirm which nodes each Segment connects,
  // in the order the line elements appear in the mesh file:
  // (1,0)-(0.75,0), (0.75,0)-(0.5,0), (0.5,0)-(0.25,0), (0.25,0)-(0,0).
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_connectivity.msh", kOneDGmsh));
  GaussLobattoLegendre gll(3);
  mesh.CreateInteriorElementNodes(gll);

  ASSERT_EQ(mesh.nodes().size(), 9u);
  const std::vector<double> expected_midpoints = {0.875, 0.625, 0.375, 0.125};
  for (size_t i = 0; i < expected_midpoints.size(); i++)
    EXPECT_NEAR(mesh.nodes().at(5 + i).x, expected_midpoints.at(i),
                EXP_NEAR_TOLERANCE)
        << "Midpoint of element " << i << " is unexpected.";
}

TEST(MeshGMSHTest, CreateInteriorElementNodesPopulatesElementNodeIDs) {
  // Element::node_ids_ is populated as a side effect of CreateInteriorNodes;
  // each of the 4 segments should end up with all 3 of its GLL(3) node IDs.
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_node_ids.msh", kOneDGmsh));
  GaussLobattoLegendre gll(3);
  mesh.CreateInteriorElementNodes(gll);

  ASSERT_EQ(mesh.n_elements(), 4u);
  for (size_t e = 0; e < mesh.n_elements(); e++)
    EXPECT_EQ(mesh.GetElement(e).node_ids().size(), 3u)
        << "Element " << e << " does not have all of its node IDs.";
}

// ---------------------------------------------------------------------------
// Mesh::Prepare
//
// Prepare() creates interior nodes, renumbers all nodes (elements share
// nodes, so a naive renumbering pass can produce ID collisions -- see
// Mesh::RenumberNodes), then checks the result. These tests exercise the
// full pipeline on the same 4-segment chain topology used above, since that
// is exactly the shared-node case renumbering has to get right.
// ---------------------------------------------------------------------------

TEST(MeshPrepareTest, DoesNotThrowAndProducesContiguousZeroBasedIDs) {
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_prepare_ids.msh", kOneDGmsh));
  GaussLobattoLegendre gll(3);

  EXPECT_NO_THROW(mesh.Prepare(gll));

  ASSERT_EQ(mesh.nodes().size(), 9u);
  for (size_t i = 0; i < mesh.nodes().size(); i++)
    EXPECT_EQ(mesh.nodes().at(i).id, i)
        << "Node at vector position " << i << " does not carry a matching ID "
                                               "(GetNode indexes by ID).";
}

TEST(MeshPrepareTest, PreservesElementConnectivityAfterRenumbering) {
  // The 4 segments form a chain: (1.0)-(0.75)-(0.5)-(0.25)-(0.0). After
  // Prepare(), each element's shared boundary must still be the *same* ID as
  // its neighbor's, and every ID must still resolve (via GetNode) to the
  // correct physical location.
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_prepare_connectivity.msh",
                          kOneDGmsh));
  GaussLobattoLegendre gll(3);
  mesh.Prepare(gll);

  ASSERT_EQ(mesh.n_elements(), 4u);
  const std::vector<double> expected_endpoint_x = {1.0, 0.75, 0.5, 0.25, 0.0};
  const std::vector<double> expected_midpoint_x = {0.875, 0.625, 0.375, 0.125};

  size_t previous_right_id = mesh.GetElement(0).node_ids().front();
  for (size_t e = 0; e < mesh.n_elements(); e++) {
    auto ids = mesh.GetElement(e).node_ids();
    ASSERT_EQ(ids.size(), 3u) << "element " << e;
    EXPECT_EQ(ids.front(), previous_right_id)
        << "element " << e << " does not connect to the previous element.";
    EXPECT_NEAR(mesh.GetNode(ids.front()).x, expected_endpoint_x.at(e),
                EXP_NEAR_TOLERANCE)
        << "element " << e << " left endpoint";
    EXPECT_NEAR(mesh.GetNode(ids.at(1)).x, expected_midpoint_x.at(e),
                EXP_NEAR_TOLERANCE)
        << "element " << e << " midpoint";
    EXPECT_NEAR(mesh.GetNode(ids.back()).x, expected_endpoint_x.at(e + 1),
                EXP_NEAR_TOLERANCE)
        << "element " << e << " right endpoint";
    previous_right_id = ids.back();
  }
}

TEST(MeshPrepareTest, HandlesSharedNodeBetweenTwoElementsWithoutIDCollisions) {
  // A minimal, hand-built (non-gmsh) case: two elements sharing a middle
  // node, which is exactly the topology that previously produced duplicate
  // IDs (see RenumberNodes).
  Mesh mesh;
  mesh.AddNodes({MakeNode(0, 0.0, 0.0, 0.0), MakeNode(1, 1.0, 0.0, 0.0),
                MakeNode(2, 2.0, 0.0, 0.0)});
  mesh.AddElement(
      std::make_unique<Segment>(std::array<size_t, 2>{0, 1}, 0, 0, mesh));
  mesh.AddElement(
      std::make_unique<Segment>(std::array<size_t, 2>{1, 2}, 0, 0, mesh));
  GaussLobattoLegendre gll(3);

  EXPECT_NO_THROW(mesh.Prepare(gll));

  ASSERT_EQ(mesh.nodes().size(), 5u);
  std::vector<size_t> ids;
  for (const auto& node : mesh.nodes()) ids.push_back(node.id);
  std::sort(ids.begin(), ids.end());
  for (size_t i = 0; i < ids.size(); i++)
    EXPECT_EQ(ids.at(i), i) << "Node IDs are not unique and contiguous.";
}

TEST(MeshGMSHTest, ThrowsOnUnsupportedSurfaceEntities) {
  constexpr char kSurfaceMesh[] = R"(
$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 1 0
$EndEntities
)";
  EXPECT_THROW(
      Mesh mesh(
          WriteTempMesh("hummingbird_mesh_test_surface.msh", kSurfaceMesh)),
      std::runtime_error);
}

TEST(MeshGMSHTest, ThrowsOnUnsupportedElementType) {
  // element_type 3 is a 4-node quad, which is not yet supported.
  constexpr char kQuadMesh[] = R"(
$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 1 0 0
1 0 0 0 0 0 0 1 4 0
$EndEntities
$Nodes
1 2 1 2
1 1 0 2
1
2
0 0 0
1 0 0
$EndNodes
$Elements
1 1 1 1
1 1 3 1
1 1 2
$EndElements
)";
  EXPECT_THROW(
      Mesh mesh(WriteTempMesh("hummingbird_mesh_test_quad.msh", kQuadMesh)),
      std::runtime_error);
}

}  // namespace hummingbird
