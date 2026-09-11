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

class BCBank : public BankBase<BC> {
 public:
  BCBank(const json& json_input);

 private:
  void Build(const json& json_input) override;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_BANKS_BC_BANK_H_
