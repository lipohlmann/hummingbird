#include "problem/problem_base.h"

namespace hummingbird {
void ProblemBase::AssembleGlobalSystem(
    const std::vector<GlobalMatrixData>& global_matrix_data) {
  auto n_vals = global_matrix_data.size();
  arma::umat locations(2, n_vals);
  std::vector<double> values(n_vals);
  for (auto i = 0; i < n_vals; i++) {
    locations(1, i) = global_matrix_data[i].row_id;
    locations(2, i) = global_matrix_data[i].col_id;
    values[i] = global_matrix_data[i].value;
  }
  arma::SpMat<double> global_system_matrix(true, locations, values, n_vals,
                                           n_vals);
}
}  // namespace hummingbird
