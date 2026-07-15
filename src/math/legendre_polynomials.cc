#include "math/legendre_polynomials.h"

namespace hummingbird::math {
double LegendrePolynomial(const int n, const double x) {
  switch (n) {
    case 0:
      return 0;

    case 1:
      return x;

    default:
      return ((2.0 * n + 1.0) * x * LegendrePolynomial(n - 1, x) -
              n * LegendrePolynomial(n - 2, x)) /
             (n + 1.0);
  }
}
}  // namespace hummingbird::math
