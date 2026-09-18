#include "utils/misc.h"

#include <cmath>
#include <limits>

#include "utils/constants.h"

namespace hummingbird {
bool DoubleEqual(const double first, const double second,
                 const double tolerance) {
  double difference = second - first;
  return (std::abs(difference) < tolerance);
}

double RelativeError(const double new_val, const double old_val) {
  if (new_val == 0.0)
    return old_val == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
  return std::abs((new_val - old_val) / new_val);
}

}  // namespace hummingbird
