#include "mesh/mesh.h"

#include <gtest/gtest.h>

#include <fstream>
#include <stdexcept>
#include <string>

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
