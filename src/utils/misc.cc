#include "utils/misc.h"

#include <cmath>

#include "utils/constants.h"

namespace hummingbird::utils {
bool DoubleEqual(const double first, const double second,
                 const double tolerance) {
  double difference = second - first;
  return (std::abs(difference) < tolerance);
}

double RelativeError(const double new_val, const double old_val) {
  return std::abs((new_val - old_val) / new_val);
}

}  // namespace hummingbird::utils
