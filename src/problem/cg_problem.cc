#include "problem/cg_problem.h"

namespace hummingbird {

CGProblem::CGProblem(const size_t n_dofs, const size_t n_ordinates)
    : ProblemBase(n_dofs, n_ordinates) {}

std::vector<GlobalMatrixData> CGProblem::AssembleGlobalMatrixData(
    const Mesh& mesh, const GaussLobattoLegendre& gll_quad,
    const MaterialBank& material_bank, const Ordinate& ordinate) {
  std::vector<GlobalMatrixData> global_matrix_data;

  for (const auto& elem : mesh.elements()) {
    auto local_mass_matrix = elem->LocalMassMatrix(gll_quad, material_bank);
    auto local_stiffness_matrix =
        elem->LocalStiffnessMatrix(gll_quad, material_bank, ordinate);
    for (auto i = 0; i < gll_quad.n_points(); i++) {
      for (auto j = 0; j < gll_quad.n_points(); j++) {
        auto node_id = elem->node_ids().at(j);
        GlobalMatrixData gmd;
        gmd.value = local_mass_matrix(i, j) + local_stiffness_matrix(i, j);
        gmd.row_id = i + node_id - 1;
        gmd.col_id = j + node_id - 1;
        global_matrix_data.push_back(std::move(gmd));
      }
    }
  }
  return global_matrix_data;
}

std::vector<GlobalForcingData> CGProblem::AssembleGlobalForcingData(
    const GaussLobattoLegendre& gll_quad, const Mesh& mesh,
    const size_t ordinate_index) {
  std::vector<GlobalForcingData> assembled_global_forcing_data;

  for (const auto& elem : mesh.elements()) {
    auto local_forcing_vector =
        elem->LocalForcingVector(gll_quad, mesh, ordinate_index);
    for (auto i = 0; i < gll_quad.n_points(); i++) {
      auto node_id = elem->node_ids().at(i);
      GlobalForcingData gfd;
      gfd.row_id = i + node_id - 1;
      gfd.value = local_forcing_vector(i);
      assembled_global_forcing_data.push_back(std::move(gfd));
    }
  }
  return assembled_global_forcing_data;
}

void CGProblem::Apply1DBCs(const Mesh& mesh, const Ordinate& ordinate,
                           const BCBank& bc_bank, const size_t ordinate_index) {
  const auto& boundary_node_ids = mesh.boundary_node_ids();

  for (const auto boundary_node_id : boundary_node_ids) {
    auto& node = mesh.GetNode(boundary_node_id);
    switch (bc_bank.GetByID(node.bc_id)) {
      case BC::VACUUM:
        double direction_dot_product =
            arma::dot(node.outward_normal, ordinate.CartesianUnitVector());
        if (direction_dot_product < 1)
          global_forcing_vectors_.at(ordinate_index)(boundary_node_id) +=
              direction_dot_product *
              solution_vectors_.at(ordinate_index)(boundary_node_id);
        else
          return;
        break;
      case BC::REFLECTIVE:
        throw std::runtime_error(
            "Reflective BCs have not been implemented for 1D yet.");
        break;

      default:
        throw std::runtime_error(
            "Invalid or unsupported BC type passed to Apply1D BCs");
        break;
    }
  }
}
}  // namespace hummingbird
