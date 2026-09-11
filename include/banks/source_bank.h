// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_PHYSICS_SOURCE_BANK_H_
#define HUMMINGBIRD_PHYSICS_SOURCE_BANK_H_

#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>

#include "physics/constant_volumetric_source.h"
#include "physics/parsed_volumetric_source.h"
#include "physics/source_base.h"

using nlohmann::json;

namespace hummingbird {

class SourceBank {
 public:
  /**
   * @brief Construct a new Source Bank object. Source structs are build using
   * the JSON input file
   *
   * @param json_input Input file as json object
   */
  SourceBank(const json& json_input);

  /**
   * @brief Get the Source object by ID
   *
   * @param id Source ID
   * @return const std::unique_ptr<SourceBase>&
   */
  const std::unique_ptr<SourceBase>& GetSource(const unsigned int id) const {
    return id_source_map_.at(id);
  }

  /**
   * @brief Get the Source object by name
   *
   * @param name Source name
   * @return const std::unique_ptr<SourceBase>&
   */
  const std::unique_ptr<SourceBase>& GetSource(const std::string& name) const {
    return GetSource(name_id_map_.at(name));
  }

 private:
  /// @brief Map from source ID to derived SourceBase object
  std::unordered_map<unsigned int, std::unique_ptr<SourceBase>> id_source_map_;

  /// @brief Map from source name to corresponding ID
  std::unordered_map<std::string, unsigned int> name_id_map_;

  /**
   * @brief Build internal maps based on JSON input
   *
   * @param json_input
   */
  void Build(const json& json_input);
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PHYSICS_SOURCE_BANK_H_
