#include "math/legendre_polynomials.h"

#include <cmath>

#include "utils/misc.h"

namespace hummingbird::math {
double LegendrePolynomial(const int n, const double x) {
  switch (n) {
    case 0:
      return 1.0;

    case 1:
      return x;

    default:
      int j = n - 1;
      return ((2.0 * static_cast<double>(j) + 1.0) * x *
                  LegendrePolynomial(j, x) -
              static_cast<double>(j) * LegendrePolynomial(j - 1, x)) /
             (static_cast<double>(j) + 1.0);
  }
}

double LegendrePolynomialDerivative(const int n, const double x) {
  if (utils::DoubleEqual(x, -1.0))
    return std::pow(-1, n - 1) * n * (n + 1) / 2.0;
  else if (utils::DoubleEqual(x, 1.0))
    return n * (n + 1) / 2.0;

  switch (n) {
    case 0:
      return 0;

    case 1:
      return 1.0;

    default:
      return (static_cast<double>(n) * LegendrePolynomial(n - 1, x) -
              static_cast<double>(n) * x * LegendrePolynomial(n, x)) /
             (1.0 - x * x);
  }
}
}  // namespace hummingbird::math
