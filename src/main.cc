#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <valarray>

#include "banks/bc_bank.h"
#include "banks/material_bank.h"
#include "banks/source_bank.h"
#include "input_parameters.h"
#include "mesh/mesh.h"
#include "problem/sem_problem.h"
#include "quadrature/angular/angular_quadrature_set.h"
#include "quadrature/gauss_lobatto_legendre.h"
#include "simulation.h"
#include "utils/json.h"
#include "utils/misc.h"
#include "utils/output.h"
#include "utils/results.h"

using nlohmann::json;
using namespace hummingbird;

int main(int argc, char** argv) {
  if (argc < 2) {
    printf("Usage: %s <input_file_name>\n", argv[0]);
    return 1;
  };

  print_header();

  std::filesystem::path working_directory(argv[1]);
  working_directory.remove_filename();

  json user_input_json = JSONFromFile(argv[1]);

  InputParams input_params = user_input_json.get<InputParams>();
  input_params.mesh_params.mesh_file =
      (working_directory / input_params.mesh_params.mesh_file).string();

  print_input_files(argv[1], input_params.mesh_params.mesh_file);

  // Create banks
  BCBank bc_bank(user_input_json);
  MaterialBank material_bank(user_input_json);
  SourceBank source_bank(user_input_json);

  // Create quadrature sets
  GaussLobattoLegendre gll_quad(input_params.sem_params.n_points);
  AngularQuadratureSet angular_quad(
      input_params.angular_treatment_params.angular_quad_set,
      input_params.angular_treatment_params.n_azim,
      input_params.angular_treatment_params.n_polar);
  size_t n_ordinates = angular_quad.get()->n_points();

  // build Mesh object
  Mesh mesh(input_params.mesh_params.mesh_file);
  mesh.Prepare(gll_quad);
  mesh.ResolveIDs(material_bank, source_bank, bc_bank);
  mesh.FindBoundaryNodes();
  mesh.SetOutwardNormals();
  mesh.InitializeNodeSolutions(n_ordinates, *angular_quad.get(), source_bank);

  // build SEMProblem
  SEMProblem sem_problem(input_params.sem_params.fe_formulation, mesh,
                         angular_quad.get()->n_points());

  Simulation simulation;

  fmt::print("Beginning source iterations.\n\n");
  print_columns();

  std::valarray<double> old_scalar_flux(0.0, sem_problem.get()->n_dofs());
  std::valarray<double> new_scalar_flux(0.0, sem_problem.get()->n_dofs());
  for (auto s_iter = 1;
       s_iter <= input_params.source_iter_params.max_iterations; s_iter++) {
    print_scatter_status(simulation.scatter_source_l2,
                         simulation.scatter_iter_error, s_iter);

    // Solve for all ordinates
    for (auto n = 0; n < n_ordinates; n++) {
      const auto ordinate = angular_quad.get()->GetAbscissa(n);

      // assemble to global system data
      auto global_matrix_data = sem_problem.get()->AssembleGlobalMatrixData(
          mesh, gll_quad, material_bank, ordinate);
      auto global_forcing_data = sem_problem.get()->AssembleGlobalForcingData(
          gll_quad, mesh, material_bank, ordinate, n);

      // form linear system
      sem_problem.get()->AssembleGlobalSystem(global_matrix_data, n);
      sem_problem.get()->AssembleGlobalForcing(global_forcing_data, n);
      sem_problem.get()->ApplyBCs(mesh, ordinate, bc_bank, n);

      // solve system
      sem_problem.get()->Solve(n);
    }
    // update scattering source
    mesh.UpdateNodeAngularFluxes(sem_problem, n_ordinates);
    mesh.UpdateNodeScalarFluxes(*angular_quad.get());
    mesh.UpdateNodeSources(*angular_quad.get(), material_bank, source_bank);

    // check source iteration convergence
    for (auto i = 0; i < mesh.n_nodes(); i++)
      new_scalar_flux[i] = mesh.GetNode(i).scalar_flux;
    std::valarray<double> error = new_scalar_flux - old_scalar_flux;

    double new_l2_error = 0;
    for (const auto& elem : mesh.elements()) {
      std::valarray<double> element_error(0.0, gll_quad.n_points());
      for (auto i = 0; i < gll_quad.n_points(); i++) {
        element_error[i] = error[elem->node_ids()[i]];
      }
      new_l2_error +=
          gll_quad.IntegrateGridFunction(element_error * element_error);
    }
    new_l2_error = std::sqrt(new_l2_error);

    simulation.scatter_iter_error =
        RelativeError(new_l2_error, simulation.scatter_source_l2);

    if (simulation.scatter_iter_error <
        input_params.source_iter_params.tolerance)
      break;

    simulation.scatter_source_l2 = new_l2_error;
    old_scalar_flux = new_scalar_flux;
    new_scalar_flux = 0.0;  // this sets all flux values to 0
  }

  // export results
  Results results(input_params.problem_params.name,
                  input_params.problem_params.output_format, mesh.nodes(),
                  *angular_quad.get());
  results.Export();

  return 0;
}
