// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_INPUT_PARAMETERS_H_
#define HUMMINGBIRD_INPUT_PARAMETERS_H_

#include <nlohmann/json.hpp>
#include <string>

#include "utils/enums.h"

using nlohmann::json;

namespace hummingbird {

struct ProblemParams {
  std::string name;
  RunMode run_mode;
  OutputFormat output_format;
};

/**
 * @brief Retrieve problem params from input
 *
 * @param j JSON object built from input file
 * @param problem_params ProblemParams struct
 */
void from_json(const json& j, ProblemParams& problem_params);

struct MeshParams {
  std::string mesh_file;
};

/**
 * @brief Retrieve mesh parameters from input
 *
 * @param j JSON object built from input file
 * @param mesh_params MeshParams struct
 */
void from_json(const json& j, MeshParams& mesh_params);

struct AngularTreatmentParams {
  AngularQuadSet angular_quad_set;
  unsigned int n_azim;
  unsigned int n_polar;
};

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

/**
 * @brief Retrieve spectral element params from input
 *
 * @param j JSON object built from input file
 * @param se_params SpectralElementParams struct
 */
void from_json(const json& j, SpectralElementParams& se_params);

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

struct InputParams {
  ProblemParams problem_params;
  MeshParams mesh_params;
  AngularTreatmentParams angular_treatment_params;
  SpectralElementParams sem_params;
  SourceIterationParams source_iter_params;
};

/**
 * @brief Retrieve input parameters from the user
 *
 * @param j JSON object built from input file
 * @param input_params InputParams struct
 */
void from_json(const json& j, InputParams& input_params);

}  // namespace hummingbird

#endif  // HUMMINGBIRD_INPUT_PARAMETERS_H_
