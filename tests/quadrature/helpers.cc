#include "../tests/quadrature/helpers.h"

namespace hummingbird::quadrature {

// Closed-form value of \int_{-1}^{1} x^k dx.
double AnalyticMonomialIntegral(const unsigned int k) {
  if (k % 2 == 1) {
    return 0.0;
  }
  return 2.0 / static_cast<double>(k + 1);
}
}  // namespace hummingbird::quadrature
