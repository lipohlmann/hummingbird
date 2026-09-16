#include "problem/problem_base.h"

#include <algorithm>
#include <format>
#include <set>
#include <stdexcept>

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

void ProblemBase::CheckGlobalMatrixData(
    const std::vector<GlobalMatrixData>& gmd) {
  // row_id/col_id pairs are expected to repeat (shared nodes get summed
  // contributions from multiple elements), so only the set of distinct IDs
  // referenced is checked for gaps.
  std::set<size_t> ids;
  for (const auto& data : gmd) {
    ids.insert(data.row_id);
    ids.insert(data.col_id);
  }

  size_t expected_id = 0;
  for (auto id : ids) {
    if (id != expected_id)
      throw std::runtime_error(std::format(
          "Global matrix data is missing node ID {}. Node IDs referenced by "
          "global matrix data must be contiguous starting from 0.",
          expected_id));
    expected_id++;
  }
}

void ProblemBase::CheckGlobalForcingData(
    const std::vector<GlobalForcingData>& gfd) {
  std::vector<size_t> ids;
  ids.reserve(gfd.size());
  for (const auto& data : gfd) ids.push_back(data.row_id);
  std::sort(ids.begin(), ids.end());

  for (size_t i = 0; i < ids.size(); i++)
    if (ids.at(i) != i)
      throw std::runtime_error(std::format(
          "Global forcing data row ID expected to be {} but was instead {}. "
          "Row IDs must be unique and contiguous starting from 0.",
          i, ids.at(i)));
}

}  // namespace hummingbird
