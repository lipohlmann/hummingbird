#include "mesh/mesh.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <format>
#include <functional>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <unordered_map>

#include "mesh/segment.h"
#include "utils/enums.h"

namespace hummingbird {

Mesh::Mesh(const std::string& msh_file) { ReadGMSH(msh_file); }

void Mesh::AddNode(const Node& node) { nodes_.push_back(node); }

void Mesh::AddNodes(const std::vector<Node>& nodes) {
  for (auto node : nodes) nodes_.push_back(std::move(node));
}

void Mesh::AddElement(std::unique_ptr<Element> element) {
  if (elements_.empty()) {
    dimension_ = element->dimension();
  } else if (element->dimension() != dimension_) {
    throw std::runtime_error(std::format(
        "Element has dimension {} but Mesh already contains elements of "
        "dimension {}. All elements in a Mesh must have the same "
        "dimension.",
        element->dimension(), dimension_));
  }
  elements_.push_back(std::move(element));
}

void Mesh::Prepare(const GaussLobattoLegendre& gll_quadrature) {
  CreateInteriorElementNodes(gll_quadrature);
  RenumberNodes();
  CheckNodeIDs();
}

void Mesh::CreateInteriorElementNodes(
    const GaussLobattoLegendre& gll_quadrature) {
  for (auto& element : elements_) {
    auto new_nodes = element->CreateInteriorNodes(nodes_, gll_quadrature);
    AddNodes(new_nodes);
  }
}

void Mesh::RenumberNodes() {
  std::unordered_map<size_t, size_t> old_to_new_id;
  std::vector<Node> renumbered_nodes(nodes_.size());

  // First, assign every old ID a new one and build the renumbered node
  // storage, without touching any element yet: elements share nodes, so an
  // ID may be seen again from a later element, but should only get a new ID
  // (and a slot in renumbered_nodes) the first time.
  size_t new_node_id = 0;
  for (auto& element : elements_) {
    for (auto old_id : element->node_ids()) {
      if (old_to_new_id.contains(old_id)) continue;
      Node node = nodes_.at(old_id);
      node.id = new_node_id;
      renumbered_nodes.at(new_node_id) = std::move(node);
      old_to_new_id.emplace(old_id, new_node_id);
      new_node_id++;
    }
  }

  // Now that every old ID maps to a final new ID, rebuild each element's
  // node ID list from scratch. Updating one ID at a time in place (e.g. via
  // SetNewNodeID) risks a new ID coincidentally colliding with a
  // not-yet-remapped old ID elsewhere in the same list.
  for (auto& element : elements_) {
    std::vector<size_t> new_ids;
    for (auto old_id : element->node_ids())
      new_ids.push_back(old_to_new_id.at(old_id));
    element->SetNodeIDs(std::move(new_ids));
  }

  nodes_ = std::move(renumbered_nodes);
}

void Mesh::CheckNodeIDs() {
  std::vector<size_t> ids;
  ids.reserve(nodes_.size());
  for (const auto& node : nodes_) ids.push_back(node.id);
  std::sort(ids.begin(), ids.end());
  for (auto i = 0; i < ids.size(); i++)
    if (ids.at(i) != i)
      throw std::runtime_error(std::format(
          "Node ID expected to be {} but was instead {}. Previous node ID: {}.",
          i, ids.at(i), ids.at(i - 1)));
}

void Mesh::ReadGMSH(const std::string& msh_file) {
  std::ifstream file(msh_file);
  if (!file.is_open())
    throw std::runtime_error(
        std::format("Could not open gmsh file: {}", msh_file));

  GmshReadState state;

  // Section headers route to their subroutine here; a new element type
  // (e.g. quads, once 2D meshes are supported) only needs a change inside
  // ReadElements, not to this dispatch.
  const std::unordered_map<std::string, std::function<void(std::ifstream&)>>
      section_handlers = {
          {"$PhysicalNames",
           [&](std::ifstream& f) { ReadPhysicalNames(f, state); }},
          {"$Entities", [&](std::ifstream& f) { ReadEntities(f, state); }},
          {"$Nodes", [&](std::ifstream& f) { ReadNodes(f, state); }},
          {"$Elements", [&](std::ifstream& f) { ReadElements(f, state); }},
      };

  std::string line;
  while (std::getline(file, line)) {
    auto handler = section_handlers.find(line);
    if (handler != section_handlers.end()) handler->second(file);
  }

  physical_names_ = std::move(state.physical_names);
}

void Mesh::ReadPhysicalNames(std::ifstream& file, GmshReadState& state) {
  size_t n_names = 0;
  file >> n_names;
  for (size_t i = 0; i < n_names; i++) {
    int dim = 0;
    int tag = 0;
    std::string name;
    file >> dim >> tag >> std::quoted(name);
    state.physical_names.emplace(tag, std::move(name));
  }
}

void Mesh::ReadEntities(std::ifstream& file, GmshReadState& state) {
  size_t n_points = 0;
  size_t n_curves = 0;
  size_t n_surfaces = 0;
  size_t n_volumes = 0;
  file >> n_points >> n_curves >> n_surfaces >> n_volumes;

  if (n_surfaces > 0 || n_volumes > 0)
    throw std::runtime_error(
        "Mesh contains surface or volume entities. Only 1D gmsh meshes "
        "(points and curves) are currently supported.");

  for (size_t i = 0; i < n_points; i++) {
    int tag = 0;
    double x = 0;
    double y = 0;
    double z = 0;
    size_t n_physical_tags = 0;
    file >> tag >> x >> y >> z >> n_physical_tags;

    std::vector<int> physical_tags(n_physical_tags);
    for (auto& physical_tag : physical_tags) file >> physical_tag;
    state.point_physical_tags.emplace(tag, std::move(physical_tags));
  }

  for (size_t i = 0; i < n_curves; i++) {
    int tag = 0;
    double min_x, min_y, min_z, max_x, max_y, max_z;
    size_t n_physical_tags = 0;
    file >> tag >> min_x >> min_y >> min_z >> max_x >> max_y >> max_z >>
        n_physical_tags;

    std::vector<int> physical_tags(n_physical_tags);
    for (auto& physical_tag : physical_tags) file >> physical_tag;
    state.curve_physical_tags.emplace(tag, std::move(physical_tags));

    size_t n_bounding_points = 0;
    file >> n_bounding_points;
    int bounding_point_tag = 0;
    for (size_t j = 0; j < n_bounding_points; j++) file >> bounding_point_tag;
  }
}

void Mesh::ReadNodes(std::ifstream& file, GmshReadState& state) {
  size_t n_blocks = 0;
  size_t n_nodes = 0;
  size_t min_tag = 0;
  size_t max_tag = 0;
  file >> n_blocks >> n_nodes >> min_tag >> max_tag;

  for (size_t block = 0; block < n_blocks; block++) {
    int entity_dim = 0;
    int entity_tag = 0;
    int parametric = 0;
    size_t n_nodes_in_block = 0;
    file >> entity_dim >> entity_tag >> parametric >> n_nodes_in_block;
    if (parametric != 0)
      throw std::runtime_error(
          "Parametric node coordinates are not currently supported.");

    std::vector<size_t> node_tags(n_nodes_in_block);
    for (auto& node_tag : node_tags) file >> node_tag;

    for (size_t i = 0; i < n_nodes_in_block; i++) {
      Node node;
      node.id = nodes_.size();
      file >> node.x >> node.y >> node.z;
      state.node_tag_to_id.emplace(node_tags.at(i), node.id);
      AddNode(node);
    }
  }
}

void Mesh::ReadElements(std::ifstream& file, GmshReadState& state) {
  // gmsh element type codes; see
  // https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
  constexpr int kPointType = 15;
  constexpr int kTwoNodeLineType = 1;

  size_t n_blocks = 0;
  size_t n_elements = 0;
  size_t min_tag = 0;
  size_t max_tag = 0;
  file >> n_blocks >> n_elements >> min_tag >> max_tag;

  for (size_t block = 0; block < n_blocks; block++) {
    int entity_dim = 0;
    int entity_tag = 0;
    int element_type = 0;
    size_t n_elements_in_block = 0;
    file >> entity_dim >> entity_tag >> element_type >> n_elements_in_block;

    if (element_type == kPointType) {
      // Point elements mark boundary entities in gmsh; they don't map to
      // their own hummingbird Element, but tag the Node they reference.
      const auto bc_id = GetBCID(state.point_physical_tags.at(entity_tag),
                                 state.physical_names);

      size_t element_tag = 0;
      size_t node_tag = 0;
      for (size_t i = 0; i < n_elements_in_block; i++) {
        file >> element_tag >> node_tag;
        Node& node = nodes_.at(state.node_tag_to_id.at(node_tag));
        node.bc_id = bc_id;
      }
      continue;
    }

    if (element_type != kTwoNodeLineType)
      throw std::runtime_error(std::format(
          "gmsh element type {} is not yet supported. Only 1D meshes "
          "(2-node line elements) can currently be read.",
          element_type));

    const auto material_id = GetMaterialID(
        state.curve_physical_tags.at(entity_tag), state.physical_names);
    const auto source_id = GetSourceID(state.curve_physical_tags.at(entity_tag),
                                       state.physical_names);

    size_t element_tag = 0;
    size_t node_tag_1 = 0;
    size_t node_tag_2 = 0;
    for (size_t i = 0; i < n_elements_in_block; i++) {
      file >> element_tag >> node_tag_1 >> node_tag_2;
      std::array<size_t, 2> boundary_node_ids = {
          state.node_tag_to_id.at(node_tag_1),
          state.node_tag_to_id.at(node_tag_2)};
      AddElement(std::make_unique<Segment>(boundary_node_ids, material_id,
                                           source_id, *this));
    }
  }
}

int Mesh::GetMaterialID(
    const std::vector<int>& curve_physical_tags,
    const std::unordered_map<int, std::string>& physical_names) const {
  for (int tag : curve_physical_tags)
    if (physical_names.at(tag).starts_with("material:")) return tag;

  throw std::runtime_error(
      "No Physical Group with a name prefixed \"material:\" was found for "
      "a curve entity.");
}

int Mesh::GetSourceID(
    const std::vector<int>& curve_physical_tags,
    const std::unordered_map<int, std::string>& physical_names) const {
  for (int tag : curve_physical_tags)
    if (physical_names.at(tag).starts_with("source:")) return tag;
  throw std::runtime_error(
      "No Physical Group with a name prefixed \"source:\" was found for "
      "a curve entity.");
}

int Mesh::GetBCID(
    const std::vector<int>& point_physical_tags,
    const std::unordered_map<int, std::string>& physical_names) const {
  for (int tag : point_physical_tags)
    if (physical_names.at(tag).starts_with("bc:")) return tag;
  throw std::runtime_error(
      "No Physical Group with a name prefixed \"bc:\" was found for "
      "a point entity.");
}

void Mesh::InitializeNodeSolutions(const size_t n_ordinates) {
  for (Node& node : nodes_) {
    node.scalar_flux = 0.0;
    node.angular_fluxes.assign(n_ordinates, 0.0);
    node.source_fluxes.assign(n_ordinates, 0.0);
  }
}

void Mesh::ResolveIDs(const MaterialBank& material_bank,
                      const SourceBank& source_bank, const BCBank& bc_bank) {
  for (auto& element : elements_) {
    const auto material_name =
        ExtractName(physical_names_.at(element->material_id()));
    try {
      element->SetMaterialID(material_bank.GetIDByName(material_name));
    } catch (const std::out_of_range&) {
      throw std::runtime_error(std::format(
          "Material \"{}\" passed to the mesh is not defined in the input "
          "file.",
          material_name));
    }

    const auto source_name =
        ExtractName(physical_names_.at(element->source_id()));
    element->SetSourceID(source_bank.GetIDByName(source_name));
  }

  for (auto& node : nodes_) {
    if (node.bc_id == 0) continue;
    const auto bc_name = ExtractName(physical_names_.at(node.bc_id));
    node.bc_id = bc_bank.GetIDByName(bc_name);
  }
}

std::string Mesh::ExtractName(const std::string& physical_name) const {
  return physical_name.substr(physical_name.find(':') + 1);
}

void Mesh::FindBoundaryNodes() {
  for (const auto& node : nodes_)
    if (node.bc_id != 0) boundary_node_ids_.push_back(node.id);
  if (dimension_ == 1) assert(boundary_node_ids_.size() == 2);
}

void Mesh::SetOutwardNormals() {
  if (dimension_ != 1)
    throw std::runtime_error(
        "SetOutwardNormals is only implemented for 1D meshes.");

  Node& node_a = nodes_.at(boundary_node_ids_.at(0));
  Node& node_b = nodes_.at(boundary_node_ids_.at(1));

  arma::vec3 direction_a = {node_a.x - node_b.x, 0.0, 0.0};
  node_a.outward_normal = arma::normalise(direction_a);
  node_b.outward_normal = -node_a.outward_normal;
}

void Mesh::UpdateNodeSources(const QuadratureBase<Ordinate>& angular_quad_set,
                             const MaterialBank material_bank,
                             const SourceBank source_bank) {
  for (auto& node : nodes_)
    UpdateSourceFluxes(node, material_bank, source_bank, angular_quad_set);
}

void Mesh::UpdateNodeScalarFluxes(
    const QuadratureBase<Ordinate>& angular_quad_set) {
  for (auto& node : nodes_) UpdateScalarFlux(node, angular_quad_set);
}

void Mesh::UpdateNodeAngularFluxes(const SEMProblem& sem_problem,
                                   const size_t n_ordinates) {
  for (auto i = 0; i < this->n_nodes(); i++) {
    for (auto n = 0; n < n_ordinates; n++) {
      nodes_[i].angular_fluxes[n] = sem_problem.get()->solution_vectors()[i][n];
    }
  }
}
}  // namespace hummingbird
