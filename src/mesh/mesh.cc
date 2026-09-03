#include "mesh/mesh.h"

#include <algorithm>
#include <format>
#include <stdexcept>

namespace hummingbird {

Mesh::Mesh(const std::string& msh_file) {}

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
}  // namespace hummingbird
