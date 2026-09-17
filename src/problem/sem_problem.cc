#include "problem/sem_problem.h"

#include <stdexcept>

namespace hummingbird {
SEMProblem::SEMProblem(const FEFormulation fe_formulation, const Mesh& mesh,
                       const size_t n_ordinates) {
  switch (fe_formulation) {
    case FEFormulation::CG:
      sem_problem = std::make_unique<CGProblem>(mesh.n_nodes(), n_ordinates);
      break;

    default:
      throw std::runtime_error("FEFormulation not supported.");
      break;
  }
}
}  // namespace hummingbird
