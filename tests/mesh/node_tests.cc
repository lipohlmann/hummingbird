// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#include "mesh/node.h"

#include <gtest/gtest.h>

#include <cmath>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "banks/material_bank.h"
#include "banks/source_bank.h"
#include "quadrature/angular/gauss_legendre_trapezoid.h"
#include "quadrature/angular/ordinate.h"
#include "utils/constants.h"

using nlohmann::json;

namespace hummingbird {

namespace {

// Helper to build a single material's JSON body.
json MakeMaterialJson(double sigma_s, double sigma_t, double sigma_f,
                      double nu) {
  return json{{"scattering_xs", sigma_s},
              {"total_xs", sigma_t},
              {"fission_xs", sigma_f},
              {"nu", nu}};
}

// Helper to build a full "materials" wrapper object, as MaterialBank expects
// to receive it.
json WrapMaterials(const std::unordered_map<std::string, json>& materials) {
  json materials_obj = json::object();
  for (const auto& [name, body] : materials) materials_obj[name] = body;
  return json{{"materials", materials_obj}};
}

// Helper to build a single "constant" source's JSON body.
json MakeConstantSourceJson(double strength) {
  return json{{"type", "constant"}, {"strength", strength}};
}

// Helper to build a single "parsed_function" source's JSON body.
json MakeParsedSourceJson(const std::string& expression) {
  return json{{"type", "parsed_function"}, {"expression", expression}};
}

// Helper to build a full "sources" wrapper object, as SourceBank expects to
// receive it.
json WrapSources(const std::unordered_map<std::string, json>& sources) {
  json sources_obj = json::object();
  for (const auto& [name, body] : sources) sources_obj[name] = body;
  return json{{"sources", sources_obj}};
}

}  // namespace

// ---------------------------------------------------------------------------
// DistanceBetweenNodes tests
// ---------------------------------------------------------------------------

TEST(DistanceBetweenNodesTest, IsZeroForCoincidentNodes) {
  Node node_1(0, 1.0, 2.0, 3.0);
  Node node_2(1, 1.0, 2.0, 3.0);
  EXPECT_NEAR(DistanceBetweenNodes(node_1, node_2), 0.0, EXP_NEAR_TOLERANCE);
}

TEST(DistanceBetweenNodesTest, ComputesEuclideanDistanceAlongSingleAxis) {
  Node node_1(0, 0.0, 0.0, 0.0);
  Node node_2(1, 5.0, 0.0, 0.0);
  EXPECT_NEAR(DistanceBetweenNodes(node_1, node_2), 5.0, EXP_NEAR_TOLERANCE);
}

TEST(DistanceBetweenNodesTest, ComputesEuclideanDistanceInThreeDimensions) {
  // A 2-3-6 right triangle in 3D has hypotenuse 7.
  Node node_1(0, 0.0, 0.0, 0.0);
  Node node_2(1, 2.0, 3.0, 6.0);
  EXPECT_NEAR(DistanceBetweenNodes(node_1, node_2), 7.0, EXP_NEAR_TOLERANCE);
}

TEST(DistanceBetweenNodesTest, IsSymmetric) {
  Node node_1(0, 1.0, -2.0, 4.0);
  Node node_2(1, -3.0, 5.0, 0.5);
  EXPECT_NEAR(DistanceBetweenNodes(node_1, node_2),
              DistanceBetweenNodes(node_2, node_1), EXP_NEAR_TOLERANCE);
}

// ---------------------------------------------------------------------------
// UpdateScalarFlux tests
// ---------------------------------------------------------------------------

TEST(UpdateScalarFluxTest, MatchesManualWeightedSumOfAngularFluxes) {
  GaussLegendreTrapezoid quad(2, 2);  // 4 quadrature points
  Node node(0, 0.0, 0.0, 0.0);
  node.angular_fluxes = {1.0, 2.0, 3.0, 4.0};

  UpdateScalarFlux(node, quad);

  double expected = 0.0;
  for (unsigned int i = 0; i < quad.n_points(); ++i)
    expected += quad.GetWeight(i) * node.angular_fluxes.at(i);
  EXPECT_NEAR(node.scalar_flux, expected, EXP_NEAR_TOLERANCE);
}

TEST(UpdateScalarFluxTest, SinglePointQuadratureScalesFluxByItsWeight) {
  GaussLegendreTrapezoid quad(1, 1);  // 1 quadrature point, weight = 4*pi
  Node node(0, 0.0, 0.0, 0.0);
  node.angular_fluxes = {2.0};

  UpdateScalarFlux(node, quad);

  EXPECT_NEAR(node.scalar_flux, quad.GetWeight(0) * 2.0, EXP_NEAR_TOLERANCE);
}

TEST(UpdateScalarFluxTest, ThrowsWhenAngularFluxesIsShorterThanQuadrature) {
  GaussLegendreTrapezoid quad(2, 2);  // 4 quadrature points
  Node node(0, 0.0, 0.0, 0.0);
  node.angular_fluxes = {1.0, 2.0};  // too few entries

  EXPECT_THROW(UpdateScalarFlux(node, quad), std::out_of_range);
}

// ---------------------------------------------------------------------------
// UpdateSourceFluxes tests
// ---------------------------------------------------------------------------

TEST(UpdateSourceFluxesTest,
    CombinesInscatteringAndConstantIndependentSourceAtEveryPoint) {
  // Material and node share id 0, since UpdateSourceFluxes looks up the
  // material bank by node.id.
  json material_input =
      WrapMaterials({{"mat", MakeMaterialJson(/*sigma_s=*/0.4,
                                              /*sigma_t=*/1.0,
                                              /*sigma_f=*/0.0, /*nu=*/0.0)}});
  MaterialBank material_bank(material_input);

  const double strength = 3.5;
  json source_input =
      WrapSources({{"src", MakeConstantSourceJson(strength)}});
  SourceBank source_bank(source_input);

  GaussLegendreTrapezoid quad(2, 2);  // 4 quadrature points
  Node node(0, 0.0, 0.0, 0.0);
  node.source_id = 0;
  node.scalar_flux = 2.0;
  node.source_fluxes.resize(quad.n_points());

  UpdateSourceFluxes(node, material_bank, source_bank, quad);

  const double expected_inscattering =
      material_bank.GetByID(0).scattering_xs / 2.0 * node.scalar_flux;
  const double expected_independent = strength / 4.0 / M_PI;
  for (unsigned int n = 0; n < quad.n_points(); ++n) {
    EXPECT_NEAR(node.source_fluxes.at(n),
                expected_inscattering + expected_independent,
                EXP_NEAR_TOLERANCE)
        << "Mismatch at quadrature point " << n;
  }
}

TEST(UpdateSourceFluxesTest, ZeroScatteringXsLeavesOnlyIndependentSource) {
  json material_input = WrapMaterials(
      {{"mat", MakeMaterialJson(/*sigma_s=*/0.0, /*sigma_t=*/1.0,
                                /*sigma_f=*/0.0, /*nu=*/0.0)}});
  MaterialBank material_bank(material_input);

  const double strength = 1.5;
  json source_input =
      WrapSources({{"src", MakeConstantSourceJson(strength)}});
  SourceBank source_bank(source_input);

  GaussLegendreTrapezoid quad(1, 1);
  Node node(0, 0.0, 0.0, 0.0);
  node.source_id = 0;
  node.scalar_flux = 10.0;  // Would matter if scattering were nonzero.
  node.source_fluxes.resize(quad.n_points());

  UpdateSourceFluxes(node, material_bank, source_bank, quad);

  EXPECT_NEAR(node.source_fluxes.at(0), strength / 4.0 / M_PI,
              EXP_NEAR_TOLERANCE);
}

TEST(UpdateSourceFluxesTest,
    UsesQuadratureAbscissaToEvaluateIndependentSourceAtEachPoint) {
  json material_input = WrapMaterials(
      {{"mat", MakeMaterialJson(/*sigma_s=*/0.0, /*sigma_t=*/1.0,
                                /*sigma_f=*/0.0, /*nu=*/0.0)}});
  MaterialBank material_bank(material_input);

  // "mu" resolves to the ordinate's x-direction cosine, so each quadrature
  // point should pick up a different independent source value.
  json source_input = WrapSources({{"src", MakeParsedSourceJson("mu")}});
  SourceBank source_bank(source_input);

  GaussLegendreTrapezoid quad(2, 2);  // 4 quadrature points
  Node node(0, 0.0, 0.0, 0.0);
  node.source_id = 0;
  node.scalar_flux = 0.0;  // Isolates the independent source contribution.
  node.source_fluxes.resize(quad.n_points());

  UpdateSourceFluxes(node, material_bank, source_bank, quad);

  for (unsigned int n = 0; n < quad.n_points(); ++n) {
    const double expected = quad.GetAbscissa(n).x();
    EXPECT_NEAR(node.source_fluxes.at(n), expected, EXP_NEAR_TOLERANCE)
        << "Mismatch at quadrature point " << n;
  }
}

}  // namespace hummingbird
