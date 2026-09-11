#include "banks/bc_bank.h"

namespace hummingbird {
BCBank::BCBank(const json& json_input) { Build(json_input); }

void BCBank::Build(const json& json_input) {
  unsigned int id = 0;
  for (const auto& [name, bc_json] :
       json_input.at("boundary_conditions").items()) {
    BC bc = bc_json.at("type");
    id_object_map_.emplace(id, std::move(bc));
    name_id_map_.emplace(name, id);
    id++;
  }
}
}  // namespace hummingbird
