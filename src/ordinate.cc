#include "ordinate.h"

#include <cmath>
#include <stdexcept>

namespace hummingbird {
Ordinate::Ordinate(const double azimuth, const double polar)
    : azimuth_(azimuth), polar_(polar) {
  CheckInput(azimuth, polar);
}

void Ordinate::CheckInput(const double azimuth, const double polar) {
  if (polar > M_PI) throw std::invalid_argument("Polar angle must be <= pi.");
  if (polar < 0) throw std::invalid_argument("Polar angle must be >= 0.");
  if (azimuth > M_PI)
    throw std::invalid_argument("Azimuthal angle must be <= pi.");
  if (azimuth < -M_PI)
    throw std::invalid_argument("Azimuthal angle must be >= 0.");
}

double Ordinate::XCosine() { return std::cos(azimuth_) * std::sin(polar_); }

double Ordinate::YCosine() { return std::sin(azimuth_) * std::sin(polar_); }

double Ordinate::ZCosine() { return std::cos(polar_); }
}  // namespace hummingbird
