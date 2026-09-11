// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_PHYSICS_SOURCE_BANK_H_
#define HUMMINGBIRD_PHYSICS_SOURCE_BANK_H_

#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>

#include "banks/bank_base.h"
#include "physics/constant_volumetric_source.h"
#include "physics/parsed_volumetric_source.h"
#include "physics/source_base.h"

using nlohmann::json;

namespace hummingbird {

class SourceBank : public BankBase<std::unique_ptr<SourceBase>> {
 public:
  /**
   * @brief Construct a new Source Bank object. Source structs are build using
   * the JSON input file
   *
   * @param json_input Input file as json object
   */
  SourceBank(const json& json_input);

 private:
  /**
   * @brief Build internal maps based on JSON input
   *
   * @param json_input
   */
  void Build(const json& json_input);
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PHYSICS_SOURCE_BANK_H_
