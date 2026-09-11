#include "input_parameters.h"

namespace hummingbird {
void from_json(const json& j, InputParams& input_params) {
  input_params.problem_params = j.get<ProblemParams>();
  input_params.mesh_params = j.get<MeshParams>();
  input_params.angular_treatment_params = j.get<AngularTreatmentParams>();
  input_params.sem_params = j.get<SpectralElementParams>();
  input_params.source_iter_params = j.get<SourceIterationParams>();
}

void from_json(const json& j, ProblemParams& problem_params) {
  j.at("problem").at("name").get_to(problem_params.name);
  j.at("problem").at("mode").get_to(problem_params.run_mode);
  j.at("problem").at("output_format").get_to(problem_params.output_format);
}

void from_json(const json& j, MeshParams& mesh_params) {
  j.at("mesh").at("filename").get_to(mesh_params.mesh_file);
}

void from_json(const json& j, AngularTreatmentParams& angular_treatment_params) {
  j.at("angular_treatment")
      .at("quadrature_set")
      .get_to(angular_treatment_params.angular_quad_set);
  j.at("angular_treatment")
      .at("n_azimuthal")
      .get_to(angular_treatment_params.n_azim);
  j.at("angular_treatment").at("n_polar").get_to(angular_treatment_params.n_polar);
}

void from_json(const json& j, SpectralElementParams& se_params) {
  j.at("spectral_elements")
      .at("transport_form")
      .get_to(se_params.transport_form);
  j.at("spectral_elements")
      .at("fe_formulation")
      .get_to(se_params.fe_formulation);
  j.at("spectral_elements").at("gll_order").get_to(se_params.gll_order);
}

void from_json(const json& j, SourceIterationParams& source_iter_params) {
  j.at("source_iteration").at("tolerance").get_to(source_iter_params.tolerance);
  j.at("source_iteration")
      .at("max_iterations")
      .get_to(source_iter_params.max_iterations);
}

}  // namespace hummingbird
