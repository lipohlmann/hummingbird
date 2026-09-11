// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#include "input_parameters.h"

#include <gtest/gtest.h>

#include <nlohmann/json.hpp>
#include <string>

#include "utils/enums.h"

using nlohmann::json;

namespace hummingbird {

// ---------------------------------------------------------------------------
// Helpers: build the JSON section(s) each from_json overload expects. Each
// returns a top-level wrapper object (e.g. {"problem": {...}}) since that's
// what the corresponding from_json calls j.at("<section>") on.
// ---------------------------------------------------------------------------

json MakeProblemJson(const std::string& name = "test_problem",
                      const std::string& mode = "fixed_source",
                      const std::string& output_format = "vtk") {
  return json{{"problem",
               {{"name", name},
                {"mode", mode},
                {"output_format", output_format}}}};
}

json MakeMeshJson(const std::string& filename = "mesh.msh") {
  return json{{"mesh", {{"filename", filename}}}};
}

json MakeAngularTreatmentJson(
    const std::string& quadrature_set = "gauss_legendre",
    unsigned int n_azimuthal = 1, unsigned int n_polar = 4) {
  return json{{"angular_treatment",
               {{"quadrature_set", quadrature_set},
                {"n_azimuthal", n_azimuthal},
                {"n_polar", n_polar}}}};
}

json MakeSpectralElementsJson(const std::string& transport_form = "saaf",
                               const std::string& fe_formulation = "cg",
                               unsigned int gll_order = 4) {
  return json{{"spectral_elements",
               {{"transport_form", transport_form},
                {"fe_formulation", fe_formulation},
                {"gll_order", gll_order}}}};
}

json MakeSourceIterationJson(double tolerance = 1e-8,
                              unsigned int max_iterations = 1000) {
  return json{{"source_iteration",
               {{"tolerance", tolerance},
                {"max_iterations", max_iterations}}}};
}

// Merges every section above into one object, as the full input file
// InputParams::from_json expects.
json MakeFullInputJson() {
  json input = json::object();
  input.update(MakeProblemJson());
  input.update(MakeMeshJson());
  input.update(MakeAngularTreatmentJson());
  input.update(MakeSpectralElementsJson());
  input.update(MakeSourceIterationJson());
  return input;
}

// ---------------------------------------------------------------------------
// ProblemParams
// ---------------------------------------------------------------------------

TEST(ProblemParamsTest, ParsesAllFields) {
  json input = MakeProblemJson("case_1", "fixed_source", "csv");

  ProblemParams params = input.get<ProblemParams>();

  EXPECT_EQ(params.name, "case_1");
  EXPECT_EQ(params.run_mode, RunMode::FIXED_SOURCE);
  EXPECT_EQ(params.output_format, OutputFormat::CSV);
}

TEST(ProblemParamsTest, ParsesVtkOutputFormat) {
  json input = MakeProblemJson("case_1", "fixed_source", "vtk");

  ProblemParams params = input.get<ProblemParams>();

  EXPECT_EQ(params.output_format, OutputFormat::VTK);
}

TEST(ProblemParamsTest, ThrowsWhenProblemKeyMissing) {
  json input = json::object();
  EXPECT_THROW(input.get<ProblemParams>(), json::out_of_range);
}

TEST(ProblemParamsTest, ThrowsWhenNameMissing) {
  json input =
      json{{"problem", {{"mode", "fixed_source"}, {"output_format", "vtk"}}}};
  EXPECT_THROW(input.get<ProblemParams>(), json::out_of_range);
}

TEST(ProblemParamsTest, ThrowsWhenModeMissing) {
  json input = json{{"problem", {{"name", "n"}, {"output_format", "vtk"}}}};
  EXPECT_THROW(input.get<ProblemParams>(), json::out_of_range);
}

TEST(ProblemParamsTest, ThrowsWhenOutputFormatMissing) {
  json input = json{{"problem", {{"name", "n"}, {"mode", "fixed_source"}}}};
  EXPECT_THROW(input.get<ProblemParams>(), json::out_of_range);
}

TEST(ProblemParamsTest, UnrecognizedOutputFormatSilentlyDefaultsToCsv) {
  // OutputFormat is deserialized via NLOHMANN_JSON_SERIALIZE_ENUM, which --
  // unlike BC's hand-written from_json (utils/enums.h) -- does not throw on
  // an unmapped string. It silently falls back to the first entry in the
  // mapping list (OutputFormat::CSV), even for a string that doesn't
  // resemble any valid value.
  json input = MakeProblemJson("n", "fixed_source", "not_a_real_format");

  ProblemParams params = input.get<ProblemParams>();

  EXPECT_EQ(params.output_format, OutputFormat::CSV);
}

// ---------------------------------------------------------------------------
// MeshParams
// ---------------------------------------------------------------------------

TEST(MeshParamsTest, ParsesFilename) {
  json input = MakeMeshJson("case_mesh.msh");

  MeshParams params = input.get<MeshParams>();

  EXPECT_EQ(params.mesh_file, "case_mesh.msh");
}

TEST(MeshParamsTest, ThrowsWhenMeshKeyMissing) {
  json input = json::object();
  EXPECT_THROW(input.get<MeshParams>(), json::out_of_range);
}

TEST(MeshParamsTest, ThrowsWhenFilenameMissing) {
  json input = json{{"mesh", json::object()}};
  EXPECT_THROW(input.get<MeshParams>(), json::out_of_range);
}

// ---------------------------------------------------------------------------
// AngularTreatmentParams
// ---------------------------------------------------------------------------

TEST(AngularTreatmentParamsTest, ParsesAllFields) {
  json input = MakeAngularTreatmentJson("gauss_legendre_chebyshev", 2, 6);

  AngularTreatmentParams params = input.get<AngularTreatmentParams>();

  EXPECT_EQ(params.angular_quad_set, AngularQuadSet::GLC);
  EXPECT_EQ(params.n_azim, 2u);
  EXPECT_EQ(params.n_polar, 6u);
}

TEST(AngularTreatmentParamsTest, ParsesEachQuadratureSet) {
  EXPECT_EQ(MakeAngularTreatmentJson("gauss_legendre")
                .get<AngularTreatmentParams>()
                .angular_quad_set,
            AngularQuadSet::GL);
  EXPECT_EQ(MakeAngularTreatmentJson("gauss_legendre_chebyshev")
                .get<AngularTreatmentParams>()
                .angular_quad_set,
            AngularQuadSet::GLC);
  EXPECT_EQ(MakeAngularTreatmentJson("gauss_legendre_trapezoid")
                .get<AngularTreatmentParams>()
                .angular_quad_set,
            AngularQuadSet::GLT);
}

TEST(AngularTreatmentParamsTest, ThrowsWhenAngularTreatmentKeyMissing) {
  json input = json::object();
  EXPECT_THROW(input.get<AngularTreatmentParams>(), json::out_of_range);
}

TEST(AngularTreatmentParamsTest, ThrowsWhenQuadratureSetMissing) {
  json input =
      json{{"angular_treatment", {{"n_azimuthal", 1}, {"n_polar", 4}}}};
  EXPECT_THROW(input.get<AngularTreatmentParams>(), json::out_of_range);
}

TEST(AngularTreatmentParamsTest, ThrowsWhenNAzimuthalMissing) {
  json input = json{
      {"angular_treatment", {{"quadrature_set", "gauss_legendre"},
                              {"n_polar", 4}}}};
  EXPECT_THROW(input.get<AngularTreatmentParams>(), json::out_of_range);
}

TEST(AngularTreatmentParamsTest, ThrowsWhenNPolarMissing) {
  json input = json{
      {"angular_treatment", {{"quadrature_set", "gauss_legendre"},
                              {"n_azimuthal", 1}}}};
  EXPECT_THROW(input.get<AngularTreatmentParams>(), json::out_of_range);
}

TEST(AngularTreatmentParamsTest,
     UnrecognizedQuadratureSetSilentlyDefaultsToGlt) {
  // Same NLOHMANN_JSON_SERIALIZE_ENUM caveat as OutputFormat above: an
  // unmapped "quadrature_set" string doesn't throw, it silently becomes
  // whichever enumerator is listed first for AngularQuadSet -- GLT.
  json input = MakeAngularTreatmentJson("not_a_real_quadrature_set", 1, 4);

  AngularTreatmentParams params = input.get<AngularTreatmentParams>();

  EXPECT_EQ(params.angular_quad_set, AngularQuadSet::GLT);
}

// ---------------------------------------------------------------------------
// SpectralElementParams
// ---------------------------------------------------------------------------

TEST(SpectralElementParamsTest, ParsesAllFields) {
  json input = MakeSpectralElementsJson("saaf", "cg", 6);

  SpectralElementParams params = input.get<SpectralElementParams>();

  EXPECT_EQ(params.transport_form, TransportForm::SAAF);
  EXPECT_EQ(params.fe_formulation, FEFormulation::CG);
  EXPECT_EQ(params.gll_order, 6u);
}

TEST(SpectralElementParamsTest, ThrowsWhenSpectralElementsKeyMissing) {
  json input = json::object();
  EXPECT_THROW(input.get<SpectralElementParams>(), json::out_of_range);
}

TEST(SpectralElementParamsTest, ThrowsWhenTransportFormMissing) {
  json input =
      json{{"spectral_elements", {{"fe_formulation", "cg"}, {"gll_order", 4}}}};
  EXPECT_THROW(input.get<SpectralElementParams>(), json::out_of_range);
}

TEST(SpectralElementParamsTest, ThrowsWhenFeFormulationMissing) {
  json input =
      json{{"spectral_elements", {{"transport_form", "saaf"}, {"gll_order", 4}}}};
  EXPECT_THROW(input.get<SpectralElementParams>(), json::out_of_range);
}

TEST(SpectralElementParamsTest, ThrowsWhenGllOrderMissing) {
  json input = json{
      {"spectral_elements", {{"transport_form", "saaf"}, {"fe_formulation", "cg"}}}};
  EXPECT_THROW(input.get<SpectralElementParams>(), json::out_of_range);
}

// ---------------------------------------------------------------------------
// SourceIterationParams
//
// BUG: from_json(const json&, SourceIterationParams&) is declared in
// input_parameters.h but has no definition anywhere in the codebase (every
// other struct declared in that header has one in input_parameters.cc).
// The tests below compile but fail to LINK with an undefined reference to
// hummingbird::from_json(nlohmann::json const&, hummingbird::
// SourceIterationParams&) -- this blocks the entire test binary from
// building, not just these cases, until a definition is added.
// ---------------------------------------------------------------------------

TEST(SourceIterationParamsTest, ParsesToleranceAndMaxIterations) {
  json input = MakeSourceIterationJson(1e-6, 500);

  SourceIterationParams params = input.get<SourceIterationParams>();

  EXPECT_DOUBLE_EQ(params.tolerance, 1e-6);
  EXPECT_EQ(params.max_iterations, 500u);
}

TEST(SourceIterationParamsTest, ThrowsWhenSourceIterationKeyMissing) {
  json input = json::object();
  EXPECT_THROW(input.get<SourceIterationParams>(), json::out_of_range);
}

TEST(SourceIterationParamsTest, ThrowsWhenToleranceMissing) {
  json input = json{{"source_iteration", {{"max_iterations", 100}}}};
  EXPECT_THROW(input.get<SourceIterationParams>(), json::out_of_range);
}

TEST(SourceIterationParamsTest, ThrowsWhenMaxIterationsMissing) {
  json input = json{{"source_iteration", {{"tolerance", 1e-5}}}};
  EXPECT_THROW(input.get<SourceIterationParams>(), json::out_of_range);
}

// ---------------------------------------------------------------------------
// InputParams (full aggregate)
// ---------------------------------------------------------------------------

TEST(InputParamsTest, ParsesAllSubParams) {
  json input = MakeFullInputJson();

  InputParams params = input.get<InputParams>();

  EXPECT_EQ(params.problem_params.name, "test_problem");
  EXPECT_EQ(params.problem_params.run_mode, RunMode::FIXED_SOURCE);
  EXPECT_EQ(params.problem_params.output_format, OutputFormat::VTK);

  EXPECT_EQ(params.mesh_params.mesh_file, "mesh.msh");

  EXPECT_EQ(params.angular_treatment_params.angular_quad_set,
            AngularQuadSet::GL);
  EXPECT_EQ(params.angular_treatment_params.n_azim, 1u);
  EXPECT_EQ(params.angular_treatment_params.n_polar, 4u);

  EXPECT_EQ(params.sem_params.transport_form, TransportForm::SAAF);
  EXPECT_EQ(params.sem_params.fe_formulation, FEFormulation::CG);
  EXPECT_EQ(params.sem_params.gll_order, 4u);

  // BUG: InputParams::from_json (src/input_parameters.cc) never assigns
  // input_params.source_iter_params -- it's simply missing from the
  // function body, unlike every other member. These two checks fail
  // because the field is left default-constructed instead of populated
  // from the "source_iteration" section above.
  EXPECT_DOUBLE_EQ(params.source_iter_params.tolerance, 1e-8);
  EXPECT_EQ(params.source_iter_params.max_iterations, 1000u);
}

TEST(InputParamsTest, ThrowsWhenProblemSectionMissing) {
  json input = MakeFullInputJson();
  input.erase("problem");
  EXPECT_THROW(input.get<InputParams>(), json::out_of_range);
}

TEST(InputParamsTest, ThrowsWhenMeshSectionMissing) {
  json input = MakeFullInputJson();
  input.erase("mesh");
  EXPECT_THROW(input.get<InputParams>(), json::out_of_range);
}

TEST(InputParamsTest, ThrowsWhenAngularTreatmentSectionMissing) {
  json input = MakeFullInputJson();
  input.erase("angular_treatment");
  EXPECT_THROW(input.get<InputParams>(), json::out_of_range);
}

TEST(InputParamsTest, ThrowsWhenSpectralElementsSectionMissing) {
  json input = MakeFullInputJson();
  input.erase("spectral_elements");
  EXPECT_THROW(input.get<InputParams>(), json::out_of_range);
}

TEST(InputParamsTest, ThrowsWhenSourceIterationSectionMissing) {
  // BUG (same root cause as ParsesAllSubParams above): every other section
  // being absent raises json::out_of_range, because InputParams::from_json
  // calls that section's own from_json via j.get<...>(). "source_iteration"
  // is never read at all, so removing it goes completely unnoticed instead
  // of throwing like its siblings do.
  json input = MakeFullInputJson();
  input.erase("source_iteration");
  EXPECT_THROW(input.get<InputParams>(), json::out_of_range);
}

}  // namespace hummingbird
