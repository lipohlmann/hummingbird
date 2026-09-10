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

struct BCParams {
  BC west;
  BC east;
  BC north = BC::NONE;
  BC south = BC::NONE;
};

NLOHMANN_JSON_SERIALIZE_ENUM(BC, {{BC::VACUUM, "vacuum"},
                                  {BC::REFLECTIVE, "reflective"}})

struct ProblemParams {
  std::string name;
  RunMode run_mode;
  SourceIterationParams source_iteration_params;
};

NLOHMANN_JSON_SERIALIZE_ENUM(RunMode,
                             {{RunMode::FIXED_SOURCE, "fixed_source"}});

struct SourceIterationParams {
  double tolerance;
  unsigned int max_iterations;
};

struct AngularTreatmentParams {
  AngularQuadSet angular_quad_set;
  unsigned int n_azim;
  unsigned int n_polar;
};

NLOHMANN_JSON_SERIALIZE_ENUM(AngularQuadSet,
                             {{AngularQuadSet::GLT, "gauss_legendre_trapezoid"},
                              {AngularQuadSet::GLC, "gauss_legendre_chebyshev"},
                              {AngularQuadSet::GL, "gauss_legendre"}});

struct SpectralElementParams {
  TransportForm transport_form;
  FEFormulation fe_formulation;
  unsigned int gll_order;
};

NLOHMANN_JSON_SERIALIZE_ENUM(TransportForm, {{TransportForm::SAAF, "SAAF"}});
NLOHMANN_JSON_SERIALIZE_ENUM(FEFormulation, {{FEFormulation::CG, "CG"}});

}  // namespace hummingbird

#endif  // HUMMINGBIRD_INPUT_PARAMETERS_H_