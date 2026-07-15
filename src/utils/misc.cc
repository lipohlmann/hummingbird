#include "utils/misc.h"

#include <cmath>

#include "utils/constants.h"

namespace hummingbird::utils {
bool DoubleEqual(const double first, const double second) {
  double difference = second - first;
  return (std::abs(difference) < TOLERANCE);
}
}  // namespace hummingbird::utils
