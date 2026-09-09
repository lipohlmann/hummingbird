// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_INPUT_PARAMETERS_H_
#define HUMMINGBIRD_INPUT_PARAMETERS_H_

#include <string>

#include "utils/enums.h"

namespace hummingbird {

struct ProblemParams {
  std::string name;
  RunMode run_mode;
  SourceIterationParams source_iteration_params;
};

struct SourceIterationParams {
  double tolerance;
  unsigned int max_iterations;
};

struct AngularTreatmentParams {
  AngularTreatmentType angular_treatment_type;
};

struct SpectralElementParams {
  TransportForm transport_form;
  FEFormulation fe_formulation;
  unsigned int gll_order;
};

}  // namespace hummingbird

#endif  // HUMMINGBIRD_INPUT_PARAMETERS_H_