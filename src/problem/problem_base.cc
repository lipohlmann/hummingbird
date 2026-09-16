#include "problem/problem_base.h"

namespace hummingbird {
void ProblemBase::AssembleGlobalSystem(
    const std::vector<GlobalMatrixData>& global_matrix_data,
    const size_t ordinate_index) {
  auto n_vals = global_matrix_data.size();
  arma::umat locations(2, n_vals);
  std::vector<double> values(n_vals);
  for (auto i = 0; i < n_vals; i++) {
    locations(0, i) = global_matrix_data[i].row_id;
    locations(1, i) = global_matrix_data[i].col_id;
    values[i] = global_matrix_data[i].value;
  }
  arma::SpMat<double> global_system_matrix(locations, values, n_vals, n_vals);
  global_system_matrices_.at(ordinate_index) = global_system_matrix;
}

void ProblemBase::Solve(const size_t ordinate_index) {
  solution_vectors_.at(ordinate_index) =
      arma::spsolve(global_system_matrices_.at(ordinate_index),
                    global_forcing_vectors_.at(ordinate_index));
}
}  // namespace hummingbird
