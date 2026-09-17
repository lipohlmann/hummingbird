// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_PHYSICS_MATERIAL_BANK_H_
#define HUMMINGBIRD_PHYSICS_MATERIAL_BANK_H_

#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>

#include "banks/bank_base.h"
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

/**
 * @brief Bank holding the materials defined in the input file
 *
 */
class MaterialBank : public BankBase<Material> {
 public:
  /**
   * @brief Construct a new Material Bank object. Material structs are built
   * using the json input
   *
   * @param json_input JSON object containing input file information
   */
  MaterialBank(const json& json_input);

 private:
  /**
   * @brief Builds the material structs and assigns them to the
   * id_material_map_ and name_id_map_ members
   *
   * @param json_input JSON object containing input file information
   */
  void Build(const json& json_input);

  /**
   * @brief Check that all materials in id_object_map_ have physically valid
   * cross sections
   *
   * @throw std::runtime_error if a material's total cross section is
   * negative, less than the sum of its scattering and fission cross
   * sections, or if its scattering, fission, or nu is negative
   */
  void CheckMaterials();
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PHYSICS_MATERIAL_BANK_H_
