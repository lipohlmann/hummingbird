#include "problem/cg_problem.h"

namespace hummingbird {

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
}  // namespace hummingbird
