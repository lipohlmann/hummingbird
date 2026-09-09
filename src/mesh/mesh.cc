#include "mesh/mesh.h"

#include <algorithm>
#include <array>
#include <format>
#include <functional>
#include <iomanip>
#include <stdexcept>
#include <unordered_map>

#include "mesh/segment.h"

namespace hummingbird {

Mesh::Mesh(const std::string& msh_file) { ReadGMSH(msh_file); }

void Mesh::AddNode(const Node& node) { nodes_.push_back(node); }

void Mesh::AddNodes(const std::vector<Node>& nodes) {
  for (auto node : nodes) nodes_.push_back(std::move(node));
}

void Mesh::AddElement(std::unique_ptr<Element> element) {
  elements_.push_back(std::move(element));
}

void Mesh::CreateInteriorElementNodes(
    const GaussLobattoLegendre& gll_quadrature) {
  for (auto& element : elements_) {
    auto new_nodes = element->CreateInteriorNodes(nodes_, gll_quadrature);
    AddNodes(new_nodes);
  }
}

void Mesh::RenumberNodes() {
  std::sort(nodes_.begin(), nodes_.end(),
            [](const Node& first, const Node& second) {
              return first.id < second.id;
            });

  size_t new_node_id = 0;
  for (auto& element : elements_) {
    auto node_ids = element->node_ids();
    for (auto id : node_ids) {
      // Elements will share nodes. Don't renumber if we've already set a new ID
      if (id < new_node_id) continue;
      nodes_.at(id).id = new_node_id;          // update global nodes
      element->SetNewNodeID(id, new_node_id);  // update element nodes
      new_node_id++;
    }
  }
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
      // Point elements only mark boundary entities in gmsh; they don't map
      // to their own hummingbird Element.
      size_t element_tag = 0;
      size_t node_tag = 0;
      for (size_t i = 0; i < n_elements_in_block; i++)
        file >> element_tag >> node_tag;
      continue;
    }

    if (element_type != kTwoNodeLineType)
      throw std::runtime_error(std::format(
          "gmsh element type {} is not yet supported. Only 1D meshes "
          "(2-node line elements) can currently be read.",
          element_type));

    const int material_id = GetMaterialID(
        state.curve_physical_tags.at(entity_tag), state.physical_names);

    size_t element_tag = 0;
    size_t node_tag_1 = 0;
    size_t node_tag_2 = 0;
    for (size_t i = 0; i < n_elements_in_block; i++) {
      file >> element_tag >> node_tag_1 >> node_tag_2;
      std::array<size_t, 2> boundary_node_ids = {
          state.node_tag_to_id.at(node_tag_1),
          state.node_tag_to_id.at(node_tag_2)};
      AddElement(std::make_unique<Segment>(boundary_node_ids, material_id));
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
}  // namespace hummingbird
