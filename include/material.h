#ifndef HUMMINGBIRD_MATERIAL_H_
#define HUMMINGBIRD_MATERIAL_H_

namespace hummingbird {

/**
 * @brief Defines a monoenergetic, constant-properties material
 *
 */
struct Material {
  /// @brief Total macroscopic cross section in 1/cm
  double total_xs;

  /// @brief Isotropic scattering macroscopic cross section in 1/cm
  double scattering_xs;

  /// @brief Fission macroscopic cross section in 1/cm
  double fission_xs;

  /// @brief Average number of neutrons produced per fission event
  double neutrons_per_fission;
};

}  // namespace hummingbird

#endif  // HUMMINGBIRD_MATERIAL_H_
