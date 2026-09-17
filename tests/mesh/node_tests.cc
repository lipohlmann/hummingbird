// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#include <gtest/gtest.h>

#include <cmath>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <vector>

#include "banks/material_bank.h"
#include "banks/source_bank.h"
#include "mesh/node.h"
#include "quadrature/angular/gauss_legendre_trapezoid.h"
#include "quadrature/angular/ordinate.h"
#include "utils/constants.h"

using nlohmann::json;

namespace hummingbird {

namespace {

// Helper to build a full "materials" wrapper object, as MaterialBank expects
// to receive it.
json WrapMaterials(double scattering_xs, double total_xs, double fission_xs,
                   double nu) {
  json body = {{"scattering_xs", scattering_xs},
               {"total_xs", total_xs},
               {"fission_xs", fission_xs},
               {"nu", nu}};
  return json{{"materials", json{{"m", body}}}};
}

// Helper to build a full "sources" wrapper object with a single constant
// source, as SourceBank expects to receive it.
json WrapConstantSource(double strength) {
  json body = {{"type", "constant"}, {"strength", strength}};
  return json{{"sources", json{{"s", body}}}};
}

// Helper to build a full "sources" wrapper object with a single
// parsed-function source, as SourceBank expects to receive it.
json WrapParsedSource(const std::string& expression) {
  json body = {{"type", "parsed_function"}, {"expression", expression}};
  return json{{"sources", json{{"s", body}}}};
}

}  // namespace

// ---------------------------------------------------------------------------
// DistanceBetweenNodes tests
// ---------------------------------------------------------------------------

TEST(DistanceBetweenNodesTest, CoincidentNodesHaveZeroDistance) {
  Node node_1(0, 1.0, 2.0, 3.0);
  Node node_2(1, 1.0, 2.0, 3.0);

  EXPECT_NEAR(DistanceBetweenNodes(node_1, node_2), 0.0, EXP_NEAR_TOLERANCE);
}

TEST(DistanceBetweenNodesTest, SingleAxisDistanceMatchesCoordinateDifference) {
  Node node_1(0, 0.0, 0.0, 0.0);
  Node node_2(1, 5.0, 0.0, 0.0);

  EXPECT_NEAR(DistanceBetweenNodes(node_1, node_2), 5.0, EXP_NEAR_TOLERANCE);
}

TEST(DistanceBetweenNodesTest, ThreeDimensionalRightTriangleDistance) {
  // A 2-3-6 right "triangle" in 3D has hypotenuse sqrt(2^2+3^2+6^2) = 7.
  Node node_1(0, 0.0, 0.0, 0.0);
  Node node_2(1, 2.0, 3.0, 6.0);

  EXPECT_NEAR(DistanceBetweenNodes(node_1, node_2), 7.0, EXP_NEAR_TOLERANCE);
}

TEST(DistanceBetweenNodesTest, DistanceIsSymmetric) {
  Node node_1(0, 1.0, -2.0, 4.0);
  Node node_2(1, -3.0, 5.0, 0.5);

  EXPECT_NEAR(DistanceBetweenNodes(node_1, node_2),
              DistanceBetweenNodes(node_2, node_1), EXP_NEAR_TOLERANCE);
}

// ---------------------------------------------------------------------------
// UpdateScalarFlux tests
// ---------------------------------------------------------------------------

TEST(UpdateScalarFluxTest, MatchesManuallyComputedQuadratureWeightedSum) {
  GaussLegendreTrapezoid quad(2, 2);

  Node node(0, 0.0, 0.0, 0.0);
  node.angular_fluxes.resize(quad.n_points());
  for (size_t i = 0; i < quad.n_points(); i++)
    node.angular_fluxes[i] = static_cast<double>(i) + 1.0;

  double expected = 0.0;
  for (size_t i = 0; i < quad.n_points(); i++)
    expected += quad.GetWeight(i) * node.angular_fluxes[i];

  UpdateScalarFlux(node, quad);

  EXPECT_NEAR(node.scalar_flux, expected, EXP_NEAR_TOLERANCE);
}

TEST(UpdateScalarFluxTest, SinglePointQuadratureScalesByFourPi) {
  // GaussLegendreTrapezoid(1, 1) has a single point whose weight is the full
  // solid angle, 4*pi.
  GaussLegendreTrapezoid quad(1, 1);

  Node node(0, 0.0, 0.0, 0.0);
  node.angular_fluxes = {2.0};

  UpdateScalarFlux(node, quad);

  EXPECT_NEAR(node.scalar_flux, 2.0 * 4.0 * M_PI, EXP_NEAR_TOLERANCE);
}

TEST(UpdateScalarFluxTest, ThrowsWhenAngularFluxesIsTooShort) {
  GaussLegendreTrapezoid quad(2, 2);

  Node node(0, 0.0, 0.0, 0.0);
  node.angular_fluxes = {1.0};  // fewer entries than quad.n_points()

  EXPECT_THROW(UpdateScalarFlux(node, quad), std::out_of_range);
}

// ---------------------------------------------------------------------------
// UpdateSourceFluxes tests
// ---------------------------------------------------------------------------

TEST(UpdateSourceFluxesTest,
     CombinesInscatteringAndConstantIndependentSourceAtEveryPoint) {
  MaterialBank material_bank(WrapMaterials(/*scattering_xs=*/1.0,
                                           /*total_xs=*/2.0,
                                           /*fission_xs=*/0.0, /*nu=*/0.0));
  SourceBank source_bank(WrapConstantSource(/*strength=*/4.0 * M_PI));
  GaussLegendreTrapezoid quad(2, 2);

  Node node(0, 0.0, 0.0, 0.0);
  node.material_id = 0;
  node.source_id = 1;  // id 0 is reserved for the zero-strength "none" source
  node.scalar_flux = 3.0;
  node.source_fluxes.resize(quad.n_points());

  UpdateSourceFluxes(node, material_bank, source_bank, quad);

  // inscattering = scattering_xs / 2 * scalar_flux = 1.0 / 2 * 3.0 = 1.5
  // independent source = strength / 4*pi = 1.0
  const double expected = 1.5 + 1.0;
  for (size_t i = 0; i < quad.n_points(); i++)
    EXPECT_NEAR(node.source_fluxes[i], expected, EXP_NEAR_TOLERANCE)
        << "Mismatch at quadrature point " << i;
}

TEST(UpdateSourceFluxesTest, ZeroScatteringIsolatesIndependentSource) {
  MaterialBank material_bank(WrapMaterials(/*scattering_xs=*/0.0,
                                           /*total_xs=*/1.0,
                                           /*fission_xs=*/0.0, /*nu=*/0.0));
  SourceBank source_bank(WrapConstantSource(/*strength=*/8.0 * M_PI));
  GaussLegendreTrapezoid quad(1, 1);

  Node node(0, 0.0, 0.0, 0.0);
  node.material_id = 0;
  node.source_id = 1;  // id 0 is reserved for the zero-strength "none" source
  node.scalar_flux = 100.0;  // should be irrelevant since scattering_xs == 0
  node.source_fluxes.resize(quad.n_points());

  UpdateSourceFluxes(node, material_bank, source_bank, quad);

  EXPECT_NEAR(node.source_fluxes[0], 2.0, EXP_NEAR_TOLERANCE);
}

TEST(UpdateSourceFluxesTest, ParsedSourcePassesEachAbscissaThrough) {
  MaterialBank material_bank(WrapMaterials(/*scattering_xs=*/0.0,
                                           /*total_xs=*/1.0,
                                           /*fission_xs=*/0.0, /*nu=*/0.0));
  // "mu" is bound to the ordinate's x-direction cosine.
  SourceBank source_bank(WrapParsedSource("mu"));
  GaussLegendreTrapezoid quad(2, 2);

  Node node(0, 0.0, 0.0, 0.0);
  node.material_id = 0;
  node.source_id = 1;  // id 0 is reserved for the zero-strength "none" source
  node.scalar_flux = 0.0;
  node.source_fluxes.resize(quad.n_points());

  UpdateSourceFluxes(node, material_bank, source_bank, quad);

  for (size_t i = 0; i < quad.n_points(); i++)
    EXPECT_NEAR(node.source_fluxes[i], quad.GetAbscissa(i).x(),
                EXP_NEAR_TOLERANCE)
        << "Mismatch at quadrature point " << i;
}

}  // namespace hummingbird
