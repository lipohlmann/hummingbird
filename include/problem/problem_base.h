// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_PROBLEM_PROBLEM_BASE_H_
#define HUMMINGBIRD_PROBLEM_PROBLEM_BASE_H_

#include <armadillo>
#include <array>

#include "mesh/mesh.h"

namespace hummingbird {
struct GlobalMatrixData {
  size_t row_id;
  size_t col_id;
  double stiffness_contrib;
  double mass_contrib;
};

class ProblemBase {
 public:
  void AssembleGlobalSystem(const GlobalMatrixData& global_matrix_data);
  virtual void Solve() = 0;

 protected:
  arma::SpMat<double> global_system_matrix_;
  arma::SpMat<double> global_mass_matrix_;
  arma::SpMat<double> global_stiffness_matrix_;
  arma::Col<double> global_forcing_vector_;
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_PROBLEM_PROBLEM_BASE_H_
