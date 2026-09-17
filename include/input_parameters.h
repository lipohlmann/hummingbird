// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_INPUT_PARAMETERS_H_
#define HUMMINGBIRD_INPUT_PARAMETERS_H_

#include <nlohmann/json.hpp>
#include <string>

#include "utils/enums.h"

using nlohmann::json;

namespace hummingbird {

/**
 * @brief Parameters from the "problem" section of the input file
 *
 */
struct ProblemParams {
  /// @brief Problem name
  std::string name;

  /// @brief Run mode for the simulation
  RunMode run_mode;

  /// @brief Format results are exported in
  OutputFormat output_format;
};

/**
 * @brief Retrieve problem params from input
 *
 * @param j JSON object built from input file
 * @param problem_params ProblemParams struct
 */
void from_json(const json& j, ProblemParams& problem_params);

/**
 * @brief Parameters from the "mesh" section of the input file
 *
 */
struct MeshParams {
  /// @brief Path to the gmsh .msh file
  std::string mesh_file;
};

/**
 * @brief Retrieve mesh parameters from input
 *
 * @param j JSON object built from input file
 * @param mesh_params MeshParams struct
 */
void from_json(const json& j, MeshParams& mesh_params);

/**
 * @brief Parameters from the "angular_treatment" section of the input file
 *
 */
struct AngularTreatmentParams {
  /// @brief Angular quadrature set to use
  AngularQuadSet angular_quad_set;

  /// @brief Number of azimuthal quadrature points
  unsigned int n_azim;

  /// @brief Number of polar quadrature points
  unsigned int n_polar;
};

/**
 * @brief Retrieve the angular treatment params from input
 *
 * @param j JSON object built from input file
 * @param angular_treatment_params AngularTreatmentParams struct
 */
void from_json(const json& j, AngularTreatmentParams& angular_treatment_params);

/**
 * @brief Parameters from the "spectral_elements" section of the input file
 *
 */
struct SpectralElementParams {
  /// @brief Form of the transport equation to solve
  TransportForm transport_form;

  /// @brief Finite element formulation to use
  FEFormulation fe_formulation;

  /// @brief Number of points per spectral element
  unsigned int n_points;
};

/**
 * @brief Retrieve spectral element params from input
 *
 * @param j JSON object built from input file
 * @param se_params SpectralElementParams struct
 */
void from_json(const json& j, SpectralElementParams& se_params);

/**
 * @brief Parameters from the "source_iteration" section of the input file
 *
 */
struct SourceIterationParams {
  /// @brief Convergence tolerance for source iteration
  double tolerance;

  /// @brief Maximum number of source iterations before giving up
  unsigned int max_iterations;
};

/**
 * @brief Retrieve the source iteration params from input
 *
 * @param j JSON object built from input file
 * @param source_iter_params SourceIterationParams struct
 */
void from_json(const json& j, SourceIterationParams& source_iter_params);

/**
 * @brief All parameters parsed from the input file
 *
 */
struct InputParams {
  /// @brief Problem parameters
  ProblemParams problem_params;

  /// @brief Mesh parameters
  MeshParams mesh_params;

  /// @brief Angular treatment parameters
  AngularTreatmentParams angular_treatment_params;

  /// @brief Spectral element parameters
  SpectralElementParams sem_params;

  /// @brief Source iteration parameters
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
