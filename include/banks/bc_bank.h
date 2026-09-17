// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_BANKS_BC_BANK_H_
#define HUMMINGBIRD_BANKS_BC_BANK_H_

#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>

#include "banks/bank_base.h"
#include "utils/enums.h"

using nlohmann::json;
namespace hummingbird {

/**
 * @brief Bank holding the boundary conditions defined in the input file
 *
 */
class BCBank : public BankBase<BC> {
 public:
  /**
   * @brief Construct a new BCBank object. BCs are built using the json input
   *
   * @param json_input JSON object containing input file information
   */
  BCBank(const json& json_input);

 private:
  /**
   * @brief Builds the BC map and assigns it to the id_object_map_ and
   * name_id_map_ members
   *
   * @param json_input JSON object containing input file information
   */
  void Build(const json& json_input) override;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_BANKS_BC_BANK_H_
