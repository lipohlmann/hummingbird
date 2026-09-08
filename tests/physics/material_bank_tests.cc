// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#include <gtest/gtest.h>

#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>

#include "physics/material.h"
#include "physics/material_bank.h"

using nlohmann::json;

namespace hummingbird {

// Helper to build a single material's JSON body.
json MakeMaterialJson(double sigma_s, double sigma_t, double sigma_f,
                      double nu) {
  return json{{"scattering_xs", sigma_s},
              {"total_xs", sigma_t},
              {"fission_xs", sigma_f},
              {"nu", nu}};
}

// Helper to build a full "materials" wrapper object, as MaterialBank expects
// to receive it (i.e. the top-level object passed to the MaterialBank
// constructor, with a "materials" key mapping names -> material bodies).
json WrapMaterials(const std::unordered_map<std::string, json>& materials) {
  json materials_obj = json::object();
  for (const auto& [name, body] : materials) materials_obj[name] = body;
  return json{{"materials", materials_obj}};
}

// ---------------------------------------------------------------------------
// from_json(Material) tests
// ---------------------------------------------------------------------------

TEST(FromJsonMaterialTest, ParsesAllFieldsCorrectly) {
  json j = {{"scattering_xs", 1.0},
            {"total_xs", 2.0},
            {"fission_xs", 0.3},
            {"nu", 2.4}};

  Material material = j.get<Material>();

  EXPECT_DOUBLE_EQ(material.scattering_xs, 1.0);
  EXPECT_DOUBLE_EQ(material.total_xs, 2.0);
  EXPECT_DOUBLE_EQ(material.fission_xs, 0.3);
  EXPECT_DOUBLE_EQ(material.nu, 2.4);
}

TEST(FromJsonMaterialTest, ThrowsWhenSigmaSMissing) {
  json j = {{"total_xs", 2.0}, {"fission_xs", 0.3}, {"nu", 2.4}};
  EXPECT_THROW(j.get<Material>(), json::out_of_range);
}

TEST(FromJsonMaterialTest, ThrowsWhenSigmaTMissing) {
  json j = {{"scattering_xs", 1.0}, {"fission_xs", 0.3}, {"nu", 2.4}};
  EXPECT_THROW(j.get<Material>(), json::out_of_range);
}

TEST(FromJsonMaterialTest, ThrowsWhenSigmaFMissing) {
  json j = {{"scattering_xs", 1.0}, {"total_xs", 2.0}, {"nu", 2.4}};
  EXPECT_THROW(j.get<Material>(), json::out_of_range);
}

TEST(FromJsonMaterialTest, ThrowsWhenNuMissing) {
  json j = {{"scattering_xs", 1.0}, {"total_xs", 2.0}, {"fission_xs", 0.3}};
  EXPECT_THROW(j.get<Material>(), json::out_of_range);
}

TEST(FromJsonMaterialTest, ThrowsOnWrongType) {
  json j = {{"scattering_xs", "not_a_number"},
            {"total_xs", 2.0},
            {"fission_xs", 0.3},
            {"nu", 2.4}};
  EXPECT_THROW(j.get<Material>(), json::type_error);
}

// ---------------------------------------------------------------------------
// MaterialBank construction / lookup tests
// ---------------------------------------------------------------------------

TEST(MaterialBankTest, ConstructsWithSingleValidMaterial) {
  json input =
      WrapMaterials({{"mms_material", MakeMaterialJson(1.0, 2.0, 0.3, 2.4)}});

  MaterialBank bank(input);

  const Material& by_id = bank.GetMaterial(0u);
  EXPECT_DOUBLE_EQ(by_id.scattering_xs, 1.0);
  EXPECT_DOUBLE_EQ(by_id.total_xs, 2.0);
  EXPECT_DOUBLE_EQ(by_id.fission_xs, 0.3);
  EXPECT_DOUBLE_EQ(by_id.nu, 2.4);

  const Material& by_name = bank.GetMaterial(std::string("mms_material"));
  EXPECT_DOUBLE_EQ(by_name.scattering_xs, 1.0);
  EXPECT_DOUBLE_EQ(by_name.total_xs, 2.0);
  EXPECT_DOUBLE_EQ(by_name.fission_xs, 0.3);
  EXPECT_DOUBLE_EQ(by_name.nu, 2.4);
}

TEST(MaterialBankTest, ConstructsWithMultipleMaterials) {
  json input = WrapMaterials({
      {"fuel", MakeMaterialJson(0.5, 1.5, 0.4, 2.5)},
      {"moderator", MakeMaterialJson(1.2, 1.2, 0.0, 0.0)},
      {"reflector", MakeMaterialJson(0.9, 1.0, 0.0, 0.0)},
  });

  MaterialBank bank(input);

  EXPECT_DOUBLE_EQ(bank.GetMaterial(0).scattering_xs, 0.5);
  EXPECT_DOUBLE_EQ(bank.GetMaterial(1).scattering_xs, 1.2);
  EXPECT_DOUBLE_EQ(bank.GetMaterial(2).scattering_xs, 0.9);

  EXPECT_DOUBLE_EQ(bank.GetMaterial(std::string("fuel")).nu, 2.5);
  EXPECT_DOUBLE_EQ(bank.GetMaterial(std::string("moderator")).total_xs, 1.2);
  EXPECT_DOUBLE_EQ(bank.GetMaterial(std::string("reflector")).total_xs, 1.0);
}

TEST(MaterialBankTest, NameAndIdMapToSameMaterialInstance) {
  json input =
      WrapMaterials({{"absorber", MakeMaterialJson(0.1, 1.0, 0.0, 0.0)}});

  MaterialBank bank(input);

  const Material& by_id = bank.GetMaterial(0);
  const Material& by_name = bank.GetMaterial(std::string("absorber"));

  EXPECT_EQ(&by_id, &by_name);
}

TEST(MaterialBankTest, GetMaterialByUnknownIdThrows) {
  json input =
      WrapMaterials({{"mms_material", MakeMaterialJson(1.0, 2.0, 0.3, 2.4)}});
  MaterialBank bank(input);

  EXPECT_THROW(bank.GetMaterial(42u), std::out_of_range);
}

TEST(MaterialBankTest, GetMaterialByUnknownNameThrows) {
  json input =
      WrapMaterials({{"mms_material", MakeMaterialJson(1.0, 2.0, 0.3, 2.4)}});
  MaterialBank bank(input);

  EXPECT_THROW(bank.GetMaterial(std::string("does_not_exist")),
               std::out_of_range);
}

TEST(MaterialBankTest, ThrowsWhenMaterialsKeyMissing) {
  json input = json::object();  // no "materials" key at all
  EXPECT_THROW(MaterialBank{input}, json::out_of_range);
}

TEST(MaterialBankTest, EmptyMaterialsObjectProducesEmptyBank) {
  json input = WrapMaterials({});
  MaterialBank bank(input);

  EXPECT_THROW(bank.GetMaterial(0u), std::out_of_range);
  EXPECT_THROW(bank.GetMaterial(std::string("anything")), std::out_of_range);
}

// ---------------------------------------------------------------------------
// CheckMaterials validation tests
// ---------------------------------------------------------------------------

TEST(MaterialBankValidationTest, AcceptsTotalXsExactlyEqualToSum) {
  // scattering(1.0) + fission(0.5) == total(1.5): should NOT throw, since
  // CheckMaterials only rejects strictly-less-than.
  json input = WrapMaterials({{"m", MakeMaterialJson(1.0, 1.5, 0.5, 1.0)}});
  EXPECT_NO_THROW(MaterialBank bank(input));
}

TEST(MaterialBankValidationTest, ThrowsWhenTotalLessThanScatteringPlusFission) {
  // scattering(1.0) + fission(0.5) = 1.5 > total(1.0)
  json input = WrapMaterials({{"m", MakeMaterialJson(1.0, 1.0, 0.5, 1.0)}});
  EXPECT_THROW(MaterialBank bank(input), std::runtime_error);
}

TEST(MaterialBankValidationTest, ThrowsOnNegativeTotalXs) {
  json input = WrapMaterials({{"m", MakeMaterialJson(0.0, -1.0, 0.0, 0.0)}});
  EXPECT_THROW(MaterialBank bank(input), std::runtime_error);
}

TEST(MaterialBankValidationTest, ThrowsOnNegativeScatteringXs) {
  json input = WrapMaterials({{"m", MakeMaterialJson(-0.1, 1.0, 0.0, 0.0)}});
  EXPECT_THROW(MaterialBank bank(input), std::runtime_error);
}

TEST(MaterialBankValidationTest, ThrowsOnNegativeFissionXs) {
  json input = WrapMaterials({{"m", MakeMaterialJson(0.0, 1.0, -0.1, 0.0)}});
  EXPECT_THROW(MaterialBank bank(input), std::runtime_error);
}

TEST(MaterialBankValidationTest, ThrowsOnNegativeNu) {
  json input = WrapMaterials({{"m", MakeMaterialJson(0.0, 1.0, 0.1, -0.5)}});
  EXPECT_THROW(MaterialBank bank(input), std::runtime_error);
}

TEST(MaterialBankValidationTest, AllZeroCrossSectionsIsValid) {
  // Vacuum-like material: everything zero should be perfectly valid.
  json input =
      WrapMaterials({{"vacuum", MakeMaterialJson(0.0, 0.0, 0.0, 0.0)}});
  EXPECT_NO_THROW(MaterialBank bank(input));
}

TEST(MaterialBankValidationTest, ValidationRunsAcrossAllMaterialsNotJustFirst) {
  // First material is valid; second is invalid. Bank construction should
  // still throw because CheckMaterials examines every material.
  json input = WrapMaterials({
      {"good", MakeMaterialJson(0.5, 1.0, 0.0, 0.0)},
      {"bad", MakeMaterialJson(2.0, 1.0, 0.0, 0.0)},  // scattering > total
  });

  EXPECT_THROW(MaterialBank bank(input), std::runtime_error);
}
}  // namespace hummingbird
