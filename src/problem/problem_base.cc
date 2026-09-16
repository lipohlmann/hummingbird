#include "problem/problem_base.h"

#include <format>
#include <set>
#include <stdexcept>

namespace hummingbird {
ProblemBase::ProblemBase(const size_t n_dofs, const size_t n_ordinates) {
  arma::SpMat<double> global_system_template(n_dofs, n_dofs);
  arma::Col<double> global_vector_template(n_dofs, arma::fill::zeros);

  global_system_matrices_.assign(n_ordinates, global_system_template);
  global_forcing_vectors_.assign(n_ordinates, global_vector_template);
  solution_vectors_.assign(n_ordinates, global_vector_template);
}

void ProblemBase::AssembleGlobalSystem(
    const std::vector<GlobalMatrixData>& global_matrix_data,
    const size_t ordinate_index) {
  CheckGlobalMatrixData(global_matrix_data);

  const auto n_dofs = global_system_matrices_.at(ordinate_index).n_rows;
  auto n_vals = global_matrix_data.size();
  arma::umat locations(2, n_vals);
  arma::Col<double> values(n_vals);
  for (auto i = 0; i < n_vals; i++) {
    locations(0, i) = global_matrix_data[i].row_id;
    locations(1, i) = global_matrix_data[i].col_id;
    values[i] = global_matrix_data[i].value;
  }
  arma::SpMat<double> global_system_matrix(locations, values, n_dofs, n_dofs);
  global_system_matrices_.at(ordinate_index) = global_system_matrix;
}

void ProblemBase::AssembleGlobalForcing(
    const std::vector<GlobalForcingData>& global_forcing_data,
    const size_t ordinate_index) {
  CheckGlobalForcingData(global_forcing_data);

  const auto n_dofs = global_forcing_vectors_.at(ordinate_index).n_elem;
  arma::Col<double> global_forcing_vector(n_dofs, arma::fill::zeros);
  for (const auto& data : global_forcing_data)
    global_forcing_vector(data.row_id) += data.value;
  global_forcing_vectors_.at(ordinate_index) = std::move(global_forcing_vector);
}

void ProblemBase::Solve(const size_t ordinate_index) {
  solution_vectors_.at(ordinate_index) =
      arma::spsolve(global_system_matrices_.at(ordinate_index),
                    global_forcing_vectors_.at(ordinate_index));
}

void ProblemBase::ApplyBCs(const Mesh& mesh, const Ordinate& ordinate,
                           const BCBank& bc_bank, const size_t ordinate_index) {
  switch (mesh.dimension()) {
    case 1:
      Apply1DBCs(mesh, ordinate, bc_bank, ordinate_index);
      break;

    default:
      throw std::runtime_error(
          "Invalid mesh ID passed in ProblemBase::ApplyBCs.");
      break;
  }
}

void ProblemBase::CheckGlobalMatrixData(
    const std::vector<GlobalMatrixData>& gmd) {
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
  std::set<size_t> ids;
  for (const auto& data : gfd) ids.insert(data.row_id);

  size_t expected_id = 0;
  for (auto id : ids) {
    if (id != expected_id)
      throw std::runtime_error(std::format(
          "Global forcing data is missing row ID {}. Row IDs referenced by "
          "global forcing data must be contiguous starting from 0.",
          expected_id));
    expected_id++;
  }
}

}  // namespace hummingbird
