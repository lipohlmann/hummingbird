#include "physics/constant_volumetric_source.h"

#include <cmath>

namespace hummingbird {
ConstantVolumetricSource::ConstantVolumetricSource(const double strength)
    : strength_(strength) {}

double ConstantVolumetricSource::EvaluateAtNode(const Node& node,
                                                const Ordinate& ordinate) {
  return strength_ / 4.0 / M_PI;
}
}  // namespace hummingbird
