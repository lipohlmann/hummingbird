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

double LegendrePolynomialDerivative(const int n, const double x) {
  double denominator = 1.0 - x * x;
  double numerator = 0;
  switch (n) {
    case 0:
      numerator = 0;

    case 1:
      numerator = 1.0;

    default:
      numerator =
          n * LegendrePolynomial(n - 1, x) - n * x * LegendrePolynomial(n, x);
  }
  return numerator / denominator;
}
}  // namespace hummingbird::math
