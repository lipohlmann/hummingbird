#include "quadrature/angular/ordinate.h"

#include <cmath>
#include <stdexcept>

namespace hummingbird {
Ordinate::Ordinate(const double azimuth, const double polar)
    : azimuth_(azimuth), polar_(polar) {
  CheckInput(azimuth, polar);
}

Ordinate Ordinate::Reflect(const Ordinate& ordinate,
                           const arma::vec3 surface_normal) {
  arma::vec3 initial_direction = ordinate.CartesianUnitVector();
  arma::vec3 new_direction =
      initial_direction -
      2.0 * arma::dot(initial_direction, surface_normal) * surface_normal;
  arma::vec3 normalized_direction = arma::normalise(new_direction);
  return Ordinate(std::atan2(normalized_direction(1), normalized_direction(0)),
                  std::acos(normalized_direction(2)));
}

void Ordinate::CheckInput(const double azimuth, const double polar) {
  if (polar > M_PI) throw std::invalid_argument("Polar angle must be <= pi.");
  if (polar < 0) throw std::invalid_argument("Polar angle must be >= 0.");
  if (azimuth > M_PI)
    throw std::invalid_argument("Azimuthal angle must be <= pi.");
  if (azimuth < -M_PI)
    throw std::invalid_argument("Azimuthal angle must be >= 0.");
}

double Ordinate::x() const { return std::cos(azimuth_) * std::sin(polar_); }

double Ordinate::y() const { return std::sin(azimuth_) * std::sin(polar_); }

double Ordinate::z() const { return std::cos(polar_); }

arma::vec3 Ordinate::CartesianUnitVector() const {
  arma::vec3 vec3;
  vec3(0) = this->x();
  vec3(1) = this->y();
  vec3(2) = this->z();
  return vec3;
}
}  // namespace hummingbird
