#include "banks/material_bank.h"

namespace hummingbird {
void from_json(const json& j, Material& material) {
  material.scattering_xs = j.at("scattering_xs").get<double>();
  material.total_xs = j.at("total_xs").get<double>();
  material.fission_xs = j.at("fission_xs").get<double>();
  material.nu = j.at("nu").get<double>();
}

MaterialBank::MaterialBank(const json& input_file_json) {
  Build(input_file_json);
  CheckMaterials();
}

void MaterialBank::Build(const json& input_file_json) {
  unsigned int id = 0;
  for (const auto& [name, material_json] :
       input_file_json.at("materials").items()) {
    Material material = material_json.get<Material>();

    id_object_map_.emplace(id, std::move(material));
    name_id_map_.emplace(name, id);
    id++;
  }
}

void MaterialBank::CheckMaterials() {
  std::vector<Material> materials;
  materials.reserve(id_object_map_.size());
  for (const auto& [id, material] : id_object_map_)
    materials.push_back(material);

  for (const auto& material : materials) {
    if (material.total_xs < (material.scattering_xs + material.fission_xs))
      throw std::runtime_error(
          "Sum of scattering and fission cross sections is greater than the "
          "total cross section. Please check over material data in input "
          "file.");
    if (material.total_xs < 0)
      throw std::runtime_error(
          "Negative total cross section passed in material. Please check over "
          "material data in input file.");
    if (material.scattering_xs < 0)
      throw std::runtime_error(
          "Negative scattering cross section passed in material. Please check "
          "over material data in input file.");
    if (material.fission_xs < 0)
      throw std::runtime_error(
          "Negative fission cross section passed in material. Please check "
          "over material data in input file.");
    if (material.nu < 0)
      throw std::runtime_error(
          "Negative fission multiplication (nu) passed in material. Please "
          "check over material data in input file.");
  }
}

}  // namespace hummingbird
