#include "input_parameters.h"

namespace hummingbird {
void from_json(const json& j, InputParams& input_params) {
  input_params.problem_params = j.get<ProblemParams>();
  input_params.mesh_params = j.get<MeshParams>();
  input_params.bc_params = j.get<BCParams>();
  input_params.angular_treatment_params = j.get<AngularTreatmentParams>();
  input_params.sem_params = j.get<SpectralElementParams>();
}

void from_json(const json& j, ProblemParams& problem_params) {
  j.at("problem").at("name").get_to(problem_params.name);
  j.at("problem").at("mode").get_to(problem_params.run_mode);
  j.at("problem").at("output_format").get_to(problem_params.output_format);
}

void from_json(const json& j, MeshParams& mesh_params) {
  j.at("mesh").at("filename").get_to(mesh_params.mesh_file);
}

void from_json(const json& j, BCParams& bc_params) {
  j.at("boundary_conditions").at("west").get_to(bc_params.west);
  j.at("boundary_conditions").at("east").get_to(bc_params.east);
}

}  // namespace hummingbird
