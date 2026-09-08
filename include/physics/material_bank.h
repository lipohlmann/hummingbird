// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_PHYSICS_MATERIAL_BANK_H_
#define HUMMINGBIRD_PHYSICS_MATERIAL_BANK_H_

#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>

#include "physics/material.h"

using nlohmann::json;

namespace hummingbird {

/**
 * @brief Get the Material information from the json input
 *
 * @param j JSON input
 * @param material Material struct to populate
 */
void from_json(const json& j, Material& material);

class MaterialBank {
 public:
  /**
   * @brief Construct a new Material Bank object. Material structs are built
   * using the json input
   *
   * @param json_input JSON object containing input file information
   */
  MaterialBank(const json& json_input);

  /**
   * @brief Get the Material object by ID
   *
   * @param id Material ID
   * @return const Material&
   */
  const Material& GetMaterial(const unsigned int id) const {
    return id_material_map_.at(id);
  }

  /**
   * @brief Get the Material object by name
   *
   * @param name Material name
   * @return const Material&
   */
  const Material& GetMaterial(const std::string& name) const {
    return GetMaterial(name_id_map_.at(name));
  }

 private:
  /// @brief Map from material ID to the Material object
  std::unordered_map<unsigned int, Material> id_material_map_;

  /// @brief Map from the material name to the corresponding ID
  std::unordered_map<std::string, unsigned int> name_id_map_;

  /**
   * @brief Builds the material structs and assigns them to the id_material_map_
   * and name_id_map_ members
   *
   * @param json_input JSON object containing input file information
   */
  void Build(const json& json_input);

  void CheckMaterials();
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PHYSICS_MATERIAL_BANK_H_
