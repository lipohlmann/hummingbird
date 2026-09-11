// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#include "banks/bc_bank.h"

#include <gtest/gtest.h>

#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "utils/enums.h"

using nlohmann::json;

namespace hummingbird {

// Helper to build a single boundary condition's JSON body.
json MakeBCJson(const std::string& type) { return json{{"type", type}}; }

// Helper to build a full "boundary_conditions" wrapper object, as BCBank
// expects to receive it (i.e. the top-level object passed to the BCBank
// constructor, with a "boundary_conditions" key mapping arbitrary names ->
// BC bodies).
json WrapBCs(const std::unordered_map<std::string, json>& bcs) {
  json bcs_obj = json::object();
  for (const auto& [name, body] : bcs) bcs_obj[name] = body;
  return json{{"boundary_conditions", bcs_obj}};
}

// ---------------------------------------------------------------------------
// Construction / arbitrary-name tests
// ---------------------------------------------------------------------------

TEST(BCBankTest, ConstructsVacuumBC) {
  json input = WrapBCs({{"left_wall", MakeBCJson("vacuum")}});

  BCBank bank(input);

  EXPECT_EQ(bank.GetByID(0u), BC::VACUUM);
  EXPECT_EQ(bank.GetByName(std::string("left_wall")), BC::VACUUM);
}

TEST(BCBankTest, ConstructsReflectiveBC) {
  json input = WrapBCs({{"symmetry_plane", MakeBCJson("reflective")}});

  BCBank bank(input);

  EXPECT_EQ(bank.GetByID(0u), BC::REFLECTIVE);
  EXPECT_EQ(bank.GetByName(std::string("symmetry_plane")), BC::REFLECTIVE);
}

TEST(BCBankTest, AcceptsArbitraryBoundaryNames) {
  // The keys under "boundary_conditions" are user-chosen labels, not tied to
  // any enum, so any string (including ones unrelated to compass directions
  // or mesh sides) must work as a lookup name.
  json input = WrapBCs({
      {"north", MakeBCJson("vacuum")},
      {"my_custom_boundary_42", MakeBCJson("reflective")},
      {"!!weird name w/ spaces!!", MakeBCJson("vacuum")},
  });

  BCBank bank(input);

  EXPECT_EQ(bank.GetByName(std::string("north")), BC::VACUUM);
  EXPECT_EQ(bank.GetByName(std::string("my_custom_boundary_42")),
            BC::REFLECTIVE);
  EXPECT_EQ(bank.GetByName(std::string("!!weird name w/ spaces!!")),
            BC::VACUUM);
}

TEST(BCBankTest, ConstructsMultipleMixedBCs) {
  json input = WrapBCs({
      {"west", MakeBCJson("vacuum")},
      {"east", MakeBCJson("reflective")},
  });

  BCBank bank(input);

  EXPECT_EQ(bank.GetByName(std::string("west")), BC::VACUUM);
  EXPECT_EQ(bank.GetByName(std::string("east")), BC::REFLECTIVE);
  EXPECT_NO_THROW(bank.GetByID(0u));
  EXPECT_NO_THROW(bank.GetByID(1u));
  EXPECT_THROW(bank.GetByID(2u), std::out_of_range);
}

TEST(BCBankTest, NameAndIdMapToSameBC) {
  json input = WrapBCs({{"only_bc", MakeBCJson("vacuum")}});
  BCBank bank(input);

  EXPECT_EQ(bank.GetByID(0u), bank.GetByName(std::string("only_bc")));
}

// ---------------------------------------------------------------------------
// Unrecognized-type-name tests
//
// BC is deserialized via a hand-written from_json (utils/enums.h) that
// throws std::invalid_argument when the JSON string doesn't match any
// mapped name, rather than silently defaulting to some enumerator. These
// tests cover BC "type" names outside {"vacuum", "reflective"} -- including
// misspellings, wrong case, and "none" (BC::NONE exists in the enum but has
// no string mapping, since it isn't a valid input-file value).
// ---------------------------------------------------------------------------

TEST(BCBankTest, UnrecognizedTypeNameThrows) {
  json input = WrapBCs({{"mystery", MakeBCJson("bogus_bc_type")}});

  EXPECT_THROW(BCBank{input}, std::invalid_argument);
}

TEST(BCBankTest, TypeNameIsCaseSensitiveAndThrows) {
  json input = WrapBCs({{"n", MakeBCJson("VACUUM")}});

  EXPECT_THROW(BCBank{input}, std::invalid_argument);
}

TEST(BCBankTest, NoneTypeNameHasNoMappingAndThrows) {
  // BC::NONE has no entry in the enum's JSON string mapping, so the string
  // "none" is unrecognized just like any other arbitrary name.
  json input = WrapBCs({{"n", MakeBCJson("none")}});

  EXPECT_THROW(BCBank{input}, std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Lookup-failure tests
// ---------------------------------------------------------------------------

TEST(BCBankTest, GetByUnknownIdThrows) {
  json input = WrapBCs({{"only_bc", MakeBCJson("vacuum")}});
  BCBank bank(input);

  EXPECT_THROW(bank.GetByID(99u), std::out_of_range);
}

TEST(BCBankTest, GetByUnknownNameThrows) {
  json input = WrapBCs({{"only_bc", MakeBCJson("vacuum")}});
  BCBank bank(input);

  EXPECT_THROW(bank.GetByName(std::string("does_not_exist")),
               std::out_of_range);
}

// ---------------------------------------------------------------------------
// Input-shape / malformed-input tests
// ---------------------------------------------------------------------------

TEST(BCBankTest, ThrowsWhenBoundaryConditionsKeyMissing) {
  json input = json::object();  // no "boundary_conditions" key at all
  EXPECT_THROW(BCBank{input}, json::out_of_range);
}

TEST(BCBankTest, EmptyBoundaryConditionsObjectProducesEmptyBank) {
  json input = WrapBCs({});
  BCBank bank(input);

  EXPECT_THROW(bank.GetByID(0u), std::out_of_range);
  EXPECT_THROW(bank.GetByName(std::string("anything")), std::out_of_range);
}

TEST(BCBankTest, ThrowsWhenTypeKeyMissing) {
  json body = json::object();  // no "type" key
  json input = WrapBCs({{"bad_bc", body}});

  EXPECT_THROW(BCBank{input}, json::out_of_range);
}
}  // namespace hummingbird
