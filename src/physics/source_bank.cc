#include "physics/source_bank.h"

namespace hummingbird {

SourceBank::SourceBank(const json& json_input) { Build(json_input); }

void SourceBank::Build(const json& json_input) {
  unsigned int id = 0;
  for (const auto& [name, source_json] : json_input.at("sources").items()) {
    if (source_json.at("type").get<std::string>() == "parsed_function") {
      id_source_map_.emplace(
          id, std::make_unique<ParsedVolumetricSource>(
                  source_json.at("expression").get<std::string>()));
      name_id_map_.emplace(name, id);
      id++;
    } else if (source_json.at("type").get<std::string>() == "constant") {
      id_source_map_.emplace(id, std::make_unique<ConstantVolumetricSource>(
                                     source_json.at("strength").get<double>()));
      name_id_map_.emplace(name, id);
      id++;
    }
  }
}
}  // namespace hummingbird
