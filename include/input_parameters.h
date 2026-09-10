// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_INPUT_PARAMETERS_H_
#define HUMMINGBIRD_INPUT_PARAMETERS_H_

#include <nlohmann/json.hpp>
#include <string>

#include "utils/enums.h"

using nlohmann::json;

namespace hummingbird {

struct InputParams {
  ProblemParams problem_params;
  std::string mesh_file;
  BCParams bc_params;
  AngularTreatmentParams angular_treatment_params;
  SpectralElementParams sem_params;
};

/**
 * @brief Retrieve input parameters from the user
 *
 * @param j JSON object built from input file
 * @param input_params InputParams struct
 */
void from_json(const json& j, InputParams& input_params);

struct BCParams {
  BC west;
  BC east;
  BC north = BC::NONE;
  BC south = BC::NONE;
};

NLOHMANN_JSON_SERIALIZE_ENUM(BC, {{BC::VACUUM, "vacuum"},
                                  {BC::REFLECTIVE, "reflective"}})

/**
 * @brief Retrieve boundary condition params from input
 *
 * @param j JSON object built from input file
 * @param bc_params BCParams struct
 */
void from_json(const json& j, BCParams& bc_params);

struct ProblemParams {
  std::string name;
  RunMode run_mode;
  SourceIterationParams source_iteration_params;
  OutputFormat output_format;
};

NLOHMANN_JSON_SERIALIZE_ENUM(RunMode,
                             {{RunMode::FIXED_SOURCE, "fixed_source"}});
NLOHMANN_JSON_SERIALIZE_ENUM(OutputFormat, {{OutputFormat::CSV, "csv"},
                                            {OutputFormat::VTK, "vtk"}});

/**
 * @brief Retrieve problem params from input
 *
 * @param j JSON object built from input file
 * @param problem_params ProblemParams struct
 */
void from_json(const json& j, ProblemParams& problem_params);

struct SourceIterationParams {
  double tolerance;
  unsigned int max_iterations;
};

/**
 * @brief Retrieve the source iteration params from input
 *
 * @param j JSON object built from input file
 * @param source_iter_params SourceIterationParams struct
 */
void from_json(const json& j, SourceIterationParams& source_iter_params);

struct AngularTreatmentParams {
  AngularQuadSet angular_quad_set;
  unsigned int n_azim;
  unsigned int n_polar;
};

NLOHMANN_JSON_SERIALIZE_ENUM(AngularQuadSet,
                             {{AngularQuadSet::GLT, "gauss_legendre_trapezoid"},
                              {AngularQuadSet::GLC, "gauss_legendre_chebyshev"},
                              {AngularQuadSet::GL, "gauss_legendre"}});

/**
 * @brief Retrieve the angular treatment params from input
 *
 * @param j JSON object built from input file
 * @param angular_treatment_params AngularTreatmentParams struct
 */
void from_json(const json& j, AngularTreatmentParams& angular_treatment_params);

struct SpectralElementParams {
  TransportForm transport_form;
  FEFormulation fe_formulation;
  unsigned int gll_order;
};

NLOHMANN_JSON_SERIALIZE_ENUM(TransportForm, {{TransportForm::SAAF, "SAAF"}});
NLOHMANN_JSON_SERIALIZE_ENUM(FEFormulation, {{FEFormulation::CG, "CG"}});

/**
 * @brief Retrieve spectral element params from input
 *
 * @param j JSON object built from input file
 * @param se_params SpectralElementParams struct
 */
void from_json(const json& j, SpectralElementParams& se_params);

}  // namespace hummingbird

#endif  // HUMMINGBIRD_INPUT_PARAMETERS_H_