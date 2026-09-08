// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#include <gtest/gtest.h>

#include <cmath>
#include <nlohmann/json.hpp>
#include <set>
#include <stdexcept>
#include <string>

#include "mesh/node.h"
#include "physics/constant_volumetric_source.h"
#include "physics/parsed_volumetric_source.h"
#include "physics/source_bank.h"
#include "physics/source_base.h"
#include "quadrature/angular/ordinate.h"

using nlohmann::json;

namespace hummingbird {

// Helper to build a single "parsed_function" source's JSON body.
json MakeParsedSourceJson(const std::string& expression) {
  return json{{"type", "parsed_function"}, {"expression", expression}};
}

// Helper to build a single "constant" source's JSON body.
json MakeConstantSourceJson(double strength) {
  return json{{"type", "constant"}, {"strength", strength}};
}

// Helper to build a full "sources" wrapper object, as SourceBank expects to
// receive it (i.e. the top-level object passed to the SourceBank
// constructor, with a "sources" key mapping names -> source bodies).
json WrapSources(const std::unordered_map<std::string, json>& sources) {
  json sources_obj = json::object();
  for (const auto& [name, body] : sources) sources_obj[name] = body;
  return json{{"sources", sources_obj}};
}

// ---------------------------------------------------------------------------
// Construction / derived-type tests
// ---------------------------------------------------------------------------

TEST(SourceBankTest, ConstructsParsedFunctionSource) {
  json input = WrapSources({{"mms_source", MakeParsedSourceJson("x^2")}});

  SourceBank bank(input);

  const std::unique_ptr<SourceBase>& by_id = bank.GetSource(0u);
  const std::unique_ptr<SourceBase>& by_name =
      bank.GetSource(std::string("mms_source"));

  ASSERT_NE(by_id.get(), nullptr);
  EXPECT_NE(dynamic_cast<ParsedVolumetricSource*>(by_id.get()), nullptr);
  EXPECT_EQ(dynamic_cast<ConstantVolumetricSource*>(by_id.get()), nullptr);
  EXPECT_EQ(by_id.get(), by_name.get());
}

TEST(SourceBankTest, ConstructsConstantSource) {
  json input = WrapSources({{"flat_source", MakeConstantSourceJson(5.0)}});

  SourceBank bank(input);

  const std::unique_ptr<SourceBase>& by_id = bank.GetSource(0u);
  const std::unique_ptr<SourceBase>& by_name =
      bank.GetSource(std::string("flat_source"));

  ASSERT_NE(by_id.get(), nullptr);
  EXPECT_NE(dynamic_cast<ConstantVolumetricSource*>(by_id.get()), nullptr);
  EXPECT_EQ(dynamic_cast<ParsedVolumetricSource*>(by_id.get()), nullptr);
  EXPECT_EQ(by_id.get(), by_name.get());
}

TEST(SourceBankTest, ConstantSourceEvaluatesToStrengthOverFourPi) {
  json input = WrapSources({{"flat_source", MakeConstantSourceJson(3.5)}});
  SourceBank bank(input);

  const std::unique_ptr<SourceBase>& source = bank.GetSource(0u);
  const double expected = 3.5 / 4.0 / M_PI;

  // Dispatch through the SourceBase* to confirm the override is actually
  // reached polymorphically (not just callable on the concrete type).
  EXPECT_DOUBLE_EQ(source->EvaluateAtNode(Node(0, 0, 0, 0), Ordinate(0, 0)),
                   expected);
  EXPECT_DOUBLE_EQ(source->EvaluateAtNode(Node(1, 2, 3, 4), Ordinate(1, 1)),
                   expected);
}

TEST(SourceBankTest, ParsedSourceEvaluateAtNodeIsReachableThroughBasePointer) {
  json input = WrapSources({{"mms_source", MakeParsedSourceJson("1")}});
  SourceBank bank(input);

  const std::unique_ptr<SourceBase>& source = bank.GetSource(0u);

  // Only checking that virtual dispatch reaches ParsedVolumetricSource's
  // override without throwing; the exact numeric result of expression
  // parsing/evaluation is exercised in detail in
  // ParsedVolumetricSourceTest (test_volumetric_sources.cpp).
  EXPECT_NO_THROW(source->EvaluateAtNode(Node(0, 0, 0, 0), Ordinate(0, 0)));
}

TEST(SourceBankTest, ParsedSourceBuiltFromJsonEvaluatesCorrectly) {
  json input = WrapSources(
      {{"mms_source", MakeParsedSourceJson("x^2 + y^2 + z^2 + mu^2")}});
  SourceBank bank(input);

  const std::unique_ptr<SourceBase>& source = bank.GetSource(0u);
  Node node(0, 2, 3, 4);
  Ordinate ordinate(0.5, 0.25);

  const double expected = node.x * node.x + node.y * node.y + node.z * node.z +
                          ordinate.x() * ordinate.x();
  EXPECT_DOUBLE_EQ(source->EvaluateAtNode(node, ordinate), expected);
}

TEST(SourceBankTest, ConstructsMultipleMixedSources) {
  json input = WrapSources({
      {"parsed_src", MakeParsedSourceJson("x^2")},
      {"const_src", MakeConstantSourceJson(7.0)},
  });

  SourceBank bank(input);

  const std::unique_ptr<SourceBase>& parsed =
      bank.GetSource(std::string("parsed_src"));
  const std::unique_ptr<SourceBase>& constant =
      bank.GetSource(std::string("const_src"));

  EXPECT_NE(dynamic_cast<ParsedVolumetricSource*>(parsed.get()), nullptr);
  EXPECT_NE(dynamic_cast<ConstantVolumetricSource*>(constant.get()), nullptr);

  // Both were assigned valid, distinct ids reachable through GetSource(id).
  EXPECT_NO_THROW(bank.GetSource(0u));
  EXPECT_NO_THROW(bank.GetSource(1u));
  // Only two sources were built, so id 2 should not exist.
  EXPECT_THROW(bank.GetSource(2u), std::out_of_range);
}

TEST(SourceBankTest, NameAndIdMapToSameSourceInstance) {
  json input = WrapSources({{"only_source", MakeConstantSourceJson(1.0)}});
  SourceBank bank(input);

  const std::unique_ptr<SourceBase>& by_id = bank.GetSource(0u);
  const std::unique_ptr<SourceBase>& by_name =
      bank.GetSource(std::string("only_source"));

  EXPECT_EQ(by_id.get(), by_name.get());
}

// ---------------------------------------------------------------------------
// Lookup-failure tests
// ---------------------------------------------------------------------------

TEST(SourceBankTest, GetSourceByUnknownIdThrows) {
  json input = WrapSources({{"only_source", MakeConstantSourceJson(1.0)}});
  SourceBank bank(input);

  EXPECT_THROW(bank.GetSource(99u), std::out_of_range);
}

TEST(SourceBankTest, GetSourceByUnknownNameThrows) {
  json input = WrapSources({{"only_source", MakeConstantSourceJson(1.0)}});
  SourceBank bank(input);

  EXPECT_THROW(bank.GetSource(std::string("does_not_exist")),
               std::out_of_range);
}

// ---------------------------------------------------------------------------
// Input-shape / malformed-input tests
// ---------------------------------------------------------------------------

TEST(SourceBankTest, ThrowsWhenSourcesKeyMissing) {
  json input = json::object();  // no "sources" key at all
  EXPECT_THROW(SourceBank{input}, json::out_of_range);
}

TEST(SourceBankTest, EmptySourcesObjectProducesEmptyBank) {
  json input = WrapSources({});
  SourceBank bank(input);

  EXPECT_THROW(bank.GetSource(0u), std::out_of_range);
  EXPECT_THROW(bank.GetSource(std::string("anything")), std::out_of_range);
}

TEST(SourceBankTest, ThrowsWhenTypeKeyMissing) {
  json body = {{"expression", "x^2"}};  // no "type" key
  json input = WrapSources({{"bad_source", body}});

  EXPECT_THROW(SourceBank{input}, json::out_of_range);
}

TEST(SourceBankTest, ThrowsWhenParsedFunctionMissingExpression) {
  json body = {{"type", "parsed_function"}};  // no "expression" key
  json input = WrapSources({{"bad_source", body}});

  EXPECT_THROW(SourceBank{input}, json::out_of_range);
}

TEST(SourceBankTest, ThrowsWhenConstantMissingStrength) {
  json body = {{"type", "constant"}};  // no "strength" key
  json input = WrapSources({{"bad_source", body}});

  EXPECT_THROW(SourceBank{input}, json::out_of_range);
}

TEST(SourceBankTest, SkipsUnrecognizedSourceTypeWithoutThrowing) {
  json input = WrapSources({
      {"good_source", MakeConstantSourceJson(2.0)},
      {"mystery_source", json{{"type", "waveform"}, {"amplitude", 1.0}}},
  });

  // Construction itself must not throw despite the unrecognized entry;
  // if it did, this line would fail with an uncaught-exception error.
  SourceBank bank(input);

  // The recognized source is present and correctly typed...
  const std::unique_ptr<SourceBase>& good =
      bank.GetSource(std::string("good_source"));
  EXPECT_NE(dynamic_cast<ConstantVolumetricSource*>(good.get()), nullptr);

  // ...while the unrecognized-type entry was never registered.
  EXPECT_THROW(bank.GetSource(std::string("mystery_source")),
               std::out_of_range);
}

TEST(SourceBankTest, IdsRemainContiguousWhenUnrecognizedTypeIsInterspersed) {
  // "mystery_source" is alphabetically between the two valid names, but
  // since Build() only increments its id counter on a successful emplace,
  // the two valid sources must still end up with contiguous ids {0, 1}
  // regardless of where the skipped entry falls in iteration order.
  json input = WrapSources({
      {"a_source", MakeConstantSourceJson(1.0)},
      {"m_mystery_source", json{{"type", "waveform"}, {"amplitude", 1.0}}},
      {"z_source", MakeParsedSourceJson("x")},
  });

  SourceBank bank(input);

  EXPECT_NO_THROW(bank.GetSource(0u));
  EXPECT_NO_THROW(bank.GetSource(1u));
  EXPECT_THROW(bank.GetSource(2u), std::out_of_range);
}
}  // namespace hummingbird