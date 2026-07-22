#include "math/legendre_polynomials.h"

#include <cmath>
#include <numeric>
#include <stdexcept>

#include "utils/constants.h"
#include "utils/misc.h"

namespace hummingbird::math {
double LegendrePolynomial(const int n, const double x) {
  switch (n) {
    case 0:
      return 1.0;

    case 1:
      return x;

    default:
      double previous = 1.0;
      double current = x;
      for (double j = 1; j < n; j++) {
        double next =
            ((2.0 * j + 1.0) * x * current - j * previous) / (j + 1.0);
        previous = current;
        current = next;
      }
      return current;
  }
}

double LegendrePolynomialPrime(const int n, const double x) {
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

double LegendrePolynomialPrimePrime(const int n, const double x) {}

std::vector<double> AllLegendreRoots(const int n) {
  unsigned int n_even =
      static_cast<unsigned int>(std::floor(static_cast<double>(n) / 2.0));
  std::vector<double> computed_roots(n);

  for (auto i = 0; i < n_even; i++) {
    computed_roots[n - 1 - i] = LegendreRoot(n, i);
    computed_roots[i] = -computed_roots[n - 1 - i];
  }
  return computed_roots;
}

double LegendreRoot(const int n, const int k) {
  if (k > n) throw std::invalid_argument("k must be less than or equal to n.");
  double x_old = std::cos((4.0 * static_cast<double>(k) + 3.0) /
                          (4.0 * static_cast<double>(n) + 2.0) * M_PI);
  double error = std::numeric_limits<double>::infinity();
  unsigned int safety = 100;
  for (auto i = 0; i < safety; i++) {
    double x_new = x_old - LegendrePolynomial(n, x_old) /
                               LegendrePolynomialPrime(n, x_old);
    error = std::abs((x_new - x_old) / x_new);
    if (error < utils::TOLERANCE) return x_new;
    x_old = x_new;
  }
  throw std::runtime_error("LegendreRoot did not converge.");
}

double LegendrePrimeRoot(const int n, const int k) {}
}  // namespace hummingbird::math
