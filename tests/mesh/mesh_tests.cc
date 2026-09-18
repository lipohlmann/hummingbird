#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>

#include "mesh/element.h"
#include "mesh/mesh.h"
#include "mesh/segment.h"
#include "quadrature/gauss_lobatto_legendre.h"
#include "utils/constants.h"

using nlohmann::json;

namespace hummingbird {

namespace {

// Mirrors cases/1D/1D_MMS_1_gmsh.msh: a unit-length 1D domain split into 4
// segments, with "bc:west"/"bc:east" boundary points and a single
// "material:mms_material" / "source:mms_source" curve region.
constexpr char kOneDGmsh[] = R"(
$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
4
0 2 "bc:west"
0 3 "bc:east"
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

std::string WriteTempMesh(const std::string& name,
                          const std::string& contents) {
  std::string path = ::testing::TempDir() + name;
  std::ofstream file(path);
  file << contents;
  return path;
}

Node MakeNode(size_t id, double x, double y, double z, unsigned int bc_id = 0) {
  Node node;
  node.id = id;
  node.x = x;
  node.y = y;
  node.z = z;
  node.bc_id = bc_id;
  return node;
}

// Minimal Element stub used only to test Mesh::AddElement's dimension
// consistency guard; the other virtuals are never exercised.
class FakeElement : public Element {
 public:
  FakeElement() : Element(0, 0) {}

  std::vector<Node> CreateInteriorNodes(const std::vector<Node>&,
                                        const GaussLobattoLegendre&) override {
    return {};
  }

  arma::Col<double> LocalForcingVector(const GaussLobattoLegendre&, const Mesh&,
                                       const size_t) override {
    return {};
  }

  arma::Mat<double> LocalStiffnessMatrix(const GaussLobattoLegendre&,
                                         const MaterialBank&,
                                         const Ordinate&) override {
    return {};
  }

  arma::SpMat<double> LocalMassMatrix(const GaussLobattoLegendre&,
                                      const MaterialBank&) override {
    return {};
  }

  unsigned int dimension() const override { return 2; }
};

// Helpers for building the banks ResolveIDs resolves kOneDGmsh's raw gmsh
// tags against.
json MakeMaterialBankJson(const std::string& name) {
  return json{{"materials",
               {{name,
                 {{"scattering_xs", 0.0},
                  {"total_xs", 1.0},
                  {"fission_xs", 0.0},
                  {"nu", 0.0}}}}}};
}

json MakeSourceBankJson(const std::string& name) {
  return json{{"sources", {{name, {{"type", "constant"}, {"strength", 1.0}}}}}};
}

json MakeBCBankJson() {
  return json{
      {"boundary_conditions",
       {{"west", {{"type", "vacuum"}}}, {"east", {{"type", "vacuum"}}}}}};
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

TEST(MeshGMSHTest, BoundaryNodesHaveMaterialIDFromMaterialPhysicalGroup) {
  // Every node in kOneDGmsh is an endpoint of the single
  // "material:mms_material" curve, so all of them should pick up its raw
  // gmsh tag (4) as material_id.
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_node_material.msh", kOneDGmsh));
  ASSERT_EQ(mesh.nodes().size(), 5u);
  for (size_t i = 0; i < mesh.nodes().size(); i++)
    EXPECT_EQ(mesh.nodes().at(i).material_id, 4)
        << "Node " << i << " has an unexpected material ID.";
}

TEST(MeshGMSHTest, BoundaryNodesHaveSourceIDFromSourcePhysicalGroup) {
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_node_source.msh", kOneDGmsh));
  ASSERT_EQ(mesh.nodes().size(), 5u);
  for (size_t i = 0; i < mesh.nodes().size(); i++)
    EXPECT_EQ(mesh.nodes().at(i).source_id, 5)
        << "Node " << i << " has an unexpected source ID.";
}

TEST(MeshGMSHTest, ConstructingFromValidMeshDoesNotThrow) {
  // Regression test: CheckMaterialIDsOnNodes runs right after ReadGMSH, so
  // every node read from a valid mesh file must already have a material ID
  // by then, not only after Prepare() creates interior nodes.
  EXPECT_NO_THROW(Mesh mesh(
      WriteTempMesh("hummingbird_mesh_test_no_throw.msh", kOneDGmsh)));
}

TEST(MeshGMSHTest, BoundaryNodesHaveBCFromBCPhysicalGroup) {
  // Pre-ResolveIDs, bc_id is still the raw gmsh Physical Group tag (2 for
  // "bc:west", 3 for "bc:east" -- see ResolveIDs tests below for the
  // resolved-to-bank-ID behavior).
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_bc.msh", kOneDGmsh));
  ASSERT_EQ(mesh.nodes().size(), 5u);
  EXPECT_EQ(mesh.nodes().at(0).bc_id, 2u);
  EXPECT_EQ(mesh.nodes().at(1).bc_id, 3u);
}

TEST(MeshGMSHTest, NonBoundaryNodesHaveNoBC) {
  Mesh mesh(
      WriteTempMesh("hummingbird_mesh_test_no_bc_interior.msh", kOneDGmsh));
  ASSERT_EQ(mesh.nodes().size(), 5u);
  for (size_t i = 2; i < mesh.nodes().size(); i++)
    EXPECT_EQ(mesh.nodes().at(i).bc_id, 0u)
        << "Node " << i << " should not have a boundary condition.";
}

// ---------------------------------------------------------------------------
// Mesh::ResolveIDs
// ---------------------------------------------------------------------------

TEST(MeshResolveIDsTest, ResolvesMaterialSourceAndBCIDsToBankAssignedIDs) {
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_resolve.msh", kOneDGmsh));
  MaterialBank material_bank(MakeMaterialBankJson("mms_material"));
  SourceBank source_bank(MakeSourceBankJson("mms_source"));
  BCBank bc_bank(MakeBCBankJson());

  mesh.ResolveIDs(material_bank, source_bank, bc_bank);

  ASSERT_EQ(mesh.n_elements(), 4u);
  for (size_t i = 0; i < mesh.n_elements(); i++) {
    EXPECT_EQ(mesh.GetElement(i).material_id(),
              static_cast<int>(material_bank.GetIDByName("mms_material")))
        << "Element " << i;
    EXPECT_EQ(mesh.GetElement(i).source_id(),
              static_cast<int>(source_bank.GetIDByName("mms_source")))
        << "Element " << i;
  }

  ASSERT_EQ(mesh.nodes().size(), 5u);
  EXPECT_EQ(mesh.nodes().at(0).bc_id, bc_bank.GetIDByName("west"));
  EXPECT_EQ(mesh.nodes().at(1).bc_id, bc_bank.GetIDByName("east"));
  for (size_t i = 2; i < mesh.nodes().size(); i++)
    EXPECT_EQ(mesh.nodes().at(i).bc_id, 0u)
        << "Node " << i << " was never tagged, should stay unresolved (0).";
}

TEST(MeshResolveIDsTest, ThrowsWhenMaterialNameNotInMaterialBank) {
  // kOneDGmsh's curves are tagged "material:mms_material", but this bank
  // only defines "some_other_material".
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_resolve_bad_material.msh",
                          kOneDGmsh));
  MaterialBank material_bank(MakeMaterialBankJson("some_other_material"));
  SourceBank source_bank(MakeSourceBankJson("mms_source"));
  BCBank bc_bank(MakeBCBankJson());

  EXPECT_THROW(mesh.ResolveIDs(material_bank, source_bank, bc_bank),
               std::runtime_error);
}

// ---------------------------------------------------------------------------
// Mesh::FindBoundaryNodes / Mesh::SetOutwardNormals
// ---------------------------------------------------------------------------

TEST(MeshBoundaryTest, FindBoundaryNodesPopulatesRealNodeIDs) {
  Mesh mesh(
      WriteTempMesh("hummingbird_mesh_test_boundary_nodes.msh", kOneDGmsh));
  GaussLobattoLegendre gll(3);
  mesh.Prepare(gll);
  MaterialBank material_bank(MakeMaterialBankJson("mms_material"));
  SourceBank source_bank(MakeSourceBankJson("mms_source"));
  BCBank bc_bank(MakeBCBankJson());
  mesh.ResolveIDs(material_bank, source_bank, bc_bank);

  mesh.FindBoundaryNodes();

  // boundary_node_ids_ must hold real node ids (usable with GetNode), not
  // the nodes' bc_id values.
  ASSERT_EQ(mesh.boundary_node_ids().size(), 2u);
  for (auto node_id : mesh.boundary_node_ids()) {
    EXPECT_NE(mesh.GetNode(node_id).bc_id, 0u);
    bool is_true_endpoint =
        std::abs(mesh.GetNode(node_id).x - 0.0) < EXP_NEAR_TOLERANCE ||
        std::abs(mesh.GetNode(node_id).x - 1.0) < EXP_NEAR_TOLERANCE;
    EXPECT_TRUE(is_true_endpoint)
        << "boundary_node_ids() entry " << node_id
        << " does not resolve to one of the mesh's true endpoints.";
  }
}

TEST(MeshBoundaryDeathTest, AssertFiresWhenBoundaryNodeCountIsNotTwoIn1D) {
  // Mesh::FindBoundaryNodes guards the 1D case with assert(), which aborts
  // the process rather than throwing a C++ exception (and is compiled out
  // entirely in an NDEBUG/release build -- this repo's debug preset, used
  // by `pixi run dev`, does not define NDEBUG). That means this needs
  // EXPECT_DEATH, not EXPECT_THROW.
  Mesh mesh;
  mesh.AddNodes(
      {MakeNode(0, 0.0, 0.0, 0.0, /*bc_id=*/1), MakeNode(1, 1.0, 0.0, 0.0)});
  mesh.AddElement(
      std::make_unique<Segment>(std::array<size_t, 2>{0, 1}, 0, 0, mesh));

  EXPECT_DEATH(mesh.FindBoundaryNodes(), "");
}

TEST(MeshBoundaryTest, SetOutwardNormalsPointsAwayFromDomain) {
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_normals.msh", kOneDGmsh));
  GaussLobattoLegendre gll(3);
  mesh.Prepare(gll);
  MaterialBank material_bank(MakeMaterialBankJson("mms_material"));
  SourceBank source_bank(MakeSourceBankJson("mms_source"));
  BCBank bc_bank(MakeBCBankJson());
  mesh.ResolveIDs(material_bank, source_bank, bc_bank);
  mesh.FindBoundaryNodes();

  mesh.SetOutwardNormals();

  for (const auto& node : mesh.nodes()) {
    arma::vec3 expected;
    if (std::abs(node.x - 0.0) < EXP_NEAR_TOLERANCE)
      expected = {-1.0, 0.0, 0.0};
    else if (std::abs(node.x - 1.0) < EXP_NEAR_TOLERANCE)
      expected = {1.0, 0.0, 0.0};
    else
      continue;  // interior node, outward_normal is unspecified

    EXPECT_NEAR(node.outward_normal(0), expected(0), EXP_NEAR_TOLERANCE)
        << "Node " << node.id;
    EXPECT_NEAR(node.outward_normal(1), expected(1), EXP_NEAR_TOLERANCE)
        << "Node " << node.id;
    EXPECT_NEAR(node.outward_normal(2), expected(2), EXP_NEAR_TOLERANCE)
        << "Node " << node.id;
  }
}

TEST(MeshBoundaryTest, SetOutwardNormalsThrowsForUnimplementedDimensions) {
  Mesh mesh;  // dimension() == 0: no elements added
  EXPECT_THROW(mesh.SetOutwardNormals(), std::runtime_error);
}

// ---------------------------------------------------------------------------
// Mesh::dimension
// ---------------------------------------------------------------------------

TEST(MeshDimensionTest, IsZeroBeforeAnyElementIsAdded) {
  Mesh mesh;
  EXPECT_EQ(mesh.dimension(), 0u);
}

TEST(MeshDimensionTest, IsOneAfterAddingSegments) {
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_dimension.msh", kOneDGmsh));
  EXPECT_EQ(mesh.dimension(), 1u);
}

TEST(MeshDimensionTest, ThrowsWhenElementDimensionsMismatch) {
  Mesh mesh;
  mesh.AddNodes({MakeNode(0, 0.0, 0.0, 0.0), MakeNode(1, 1.0, 0.0, 0.0)});
  mesh.AddElement(
      std::make_unique<Segment>(std::array<size_t, 2>{0, 1}, 0, 0, mesh));
  EXPECT_THROW(mesh.AddElement(std::make_unique<FakeElement>()),
               std::runtime_error);
}

TEST(MeshGMSHTest, ThrowsWhenPointEntityHasNoBCPhysicalGroup) {
  // The point's only Physical Group name doesn't start with "bc:", so
  // Mesh::GetBCID has nothing to find.
  constexpr char kNoBCPointMesh[] = R"(
$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
1
0 2 "not_a_bc"
$EndPhysicalNames
$Entities
1 0 0 0
1 0 0 0 1 2
$EndEntities
$Nodes
1 1 1 1
0 1 0 1
1
0 0 0
$EndNodes
$Elements
1 1 1 1
0 1 15 1
1 1
$EndElements
)";
  EXPECT_THROW(Mesh mesh(WriteTempMesh("hummingbird_mesh_test_no_bc_group.msh",
                                       kNoBCPointMesh)),
               std::runtime_error);
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
    EXPECT_EQ(mesh.nodes().at(i).id, i) << "Node at vector position " << i
                                        << " does not carry a matching ID "
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

TEST(MeshPrepareTest,
     BoundaryConditionsSurviveRenumberingAndInteriorNodeCreation) {
  // The chain's two true endpoints are tagged by kOneDGmsh's point entities:
  // x=0.0 by "bc:west" (raw tag 2), x=1.0 by "bc:east" (raw tag 3). After
  // Prepare() (interior node creation + renumbering), those two nodes must
  // still carry their original bc_id, while every other node -- the
  // pre-existing curve nodes and the newly created GLL interior nodes --
  // must remain untagged (bc_id == 0).
  Mesh mesh(WriteTempMesh("hummingbird_mesh_test_prepare_bc.msh", kOneDGmsh));
  GaussLobattoLegendre gll(3);
  mesh.Prepare(gll);

  ASSERT_EQ(mesh.nodes().size(), 9u);
  for (const auto& node : mesh.nodes()) {
    unsigned int expected_bc_id = 0u;
    if (std::abs(node.x - 0.0) < EXP_NEAR_TOLERANCE) expected_bc_id = 2u;
    if (std::abs(node.x - 1.0) < EXP_NEAR_TOLERANCE) expected_bc_id = 3u;
    EXPECT_EQ(node.bc_id, expected_bc_id)
        << "Node " << node.id << " at x=" << node.x
        << " has an unexpected boundary condition.";
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
  EXPECT_THROW(Mesh mesh(WriteTempMesh("hummingbird_mesh_test_surface.msh",
                                       kSurfaceMesh)),
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
