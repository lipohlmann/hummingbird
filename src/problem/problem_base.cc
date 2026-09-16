#include "problem/problem_base.h"

namespace hummingbird {
void ProblemBase::AssembleGlobalSystem(
    const std::vector<GlobalMatrixData>& global_matrix_data,
    const size_t ordinate_index) {
  CheckGlobalMatrixData(global_matrix_data);

  auto n_vals = global_matrix_data.size();
  arma::umat locations(2, n_vals);
  arma::Col<double> values(n_vals);
  for (auto i = 0; i < n_vals; i++) {
    locations(0, i) = global_matrix_data[i].row_id;
    locations(1, i) = global_matrix_data[i].col_id;
    values[i] = global_matrix_data[i].value;
  }
  arma::SpMat<double> global_system_matrix(locations, values, n_vals, n_vals);
  global_system_matrices_.at(ordinate_index) = global_system_matrix;
}

void ProblemBase::AssembleGlobalForcing(
    const std::vector<GlobalForcingData>& global_forcing_data,
    const size_t ordinate_index) {
  CheckGlobalForcingData(global_forcing_data);

  arma::Col<double> global_forcing_vector(global_forcing_data.size(),
                                          arma::fill::zeros);
  for (const auto& data : global_forcing_data)
    global_forcing_vector(data.row_id) = data.value;
  global_forcing_vectors_.at(ordinate_index) = std::move(global_forcing_vector);
}

void ProblemBase::Solve(const size_t ordinate_index) {
  solution_vectors_.at(ordinate_index) =
      arma::spsolve(global_system_matrices_.at(ordinate_index),
                    global_forcing_vectors_.at(ordinate_index));
}

}  // namespace hummingbird
