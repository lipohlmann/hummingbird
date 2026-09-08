#include "physics/source_bank.h"

namespace hummingbird {

void SourceBank::Build(const json& json_input) {
  unsigned int id = 0;
  for (const auto& [name, source_json] : json_input.at("sources").items()) {
    if (source_json.at("type").get<std::string>() == "parsed_function") {
      ParsedVolumetricSource pvs(
          source_json.at("expression").get<std::string>());
      id_source_map_.emplace(id, std::move(pvs));
      name_id_map_.emplace(name, id);
      id++;
    } else if (source_json.at("type").get<std::string>() == "constant") {
      ConstantVolumetricSource cvs(source_json.at("strength").get<double>());
      id_source_map_.emplace(id, std::move(cvs));
      name_id_map_.emplace(name, id);
      id++;
    }
  }
}
}  // namespace hummingbird
