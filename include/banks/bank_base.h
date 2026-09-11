// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_BANKS_BANK_BASE_H_
#define HUMMINGBIRD_BANKS_BANK_BASE_H_

#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>

using nlohmann::json;

namespace hummingbird {
template <typename T>
class BankBase {
 public:
  /**
   * @brief Get the object by its ID
   *
   * @param id ID
   * @return const T&
   */
  const T& GetByID(const unsigned int id) const {
    return id_object_map_.at(id);
  }

  /**
   * @brief Get the object by its name in the JSON input file
   *
   * @param name Object name
   * @return const T&
   */
  const T& GetByName(const std::string& name) const {
    return GetByID(name_id_map_.at(name));
  }

 protected:
  /// @brief ID to object map
  std::unordered_map<unsigned int, T> id_object_map_;

  /// @brief name to object ID map
  std::unordered_map<std::string, unsigned int> name_id_map_;

  /**
   * @brief Build the bank
   *
   * @param json_input JSON input
   */
  virtual void Build(const json& json_input) = 0;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_BANKS_BANK_BASE_H_
