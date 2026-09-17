#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

#include "banks/bc_bank.h"
#include "banks/material_bank.h"
#include "banks/source_bank.h"
#include "input_parameters.h"
#include "mesh/mesh.h"
#include "quadrature/angular/angular_quadrature_set.h"
#include "quadrature/gauss_lobatto_legendre.h"
#include "simulation.h"
#include "utils/json.h"
#include "utils/output.h"

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

  // build Mesh object
  Mesh mesh(input_params.mesh_params.mesh_file);
  mesh.Prepare(gll_quad);
  mesh.ResolveIDs(material_bank, source_bank, bc_bank);
  mesh.FindBoundaryNodes();
  mesh.SetOutwardNormals();

  Simulation simulation;

  // form local matrices

  // assemble to global system

  // solve system (should be a one-liner?)

  // check source iteration convergence

  // export results

  return 0;
}
