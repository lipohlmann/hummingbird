#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "math/legendre_polynomials.h"

namespace hummingbird::math {

// Representative x values in the open interval, plus endpoints.
const std::vector<double> kInteriorX = {
    -0.9, -0.5, -0.1, 0.0, 0.1, 0.5, 1.0 / std::sqrt(2.0), 0.9};
const std::vector<double> kAllX = {
    -1.0, -0.9, -0.5, -0.1, 0.0, 0.1, 0.5, 1.0 / std::sqrt(2.0), 0.9, 1.0};

// =======================================================================
// LegendrePolynomial(n, x)
// =======================================================================

TEST(LegendrePolynomial, P0IsOne) {
  for (double x : kAllX) {
    EXPECT_DOUBLE_EQ(LegendrePolynomial(0, x), 1.0);
  }
}

TEST(LegendrePolynomial, P1IsX) {
  for (double x : kAllX) {
    EXPECT_DOUBLE_EQ(LegendrePolynomial(1, x), x);
  }
}

TEST(LegendrePolynomial, P2ClosedForm) {
  for (double x : kAllX) {
    double expected = 0.5 * (3.0 * x * x - 1.0);
    EXPECT_DOUBLE_EQ(LegendrePolynomial(2, x), expected);
  }
}

TEST(LegendrePolynomial, P3ClosedForm) {
  for (double x : kAllX) {
    double expected = 0.5 * (5.0 * x * x * x - 3.0 * x);
    EXPECT_DOUBLE_EQ(LegendrePolynomial(3, x), expected);
  }
}

TEST(LegendrePolynomial, P4ClosedForm) {
  for (double x : kAllX) {
    double x2 = x * x;
    double expected = (35.0 * x2 * x2 - 30.0 * x2 + 3.0) / 8.0;
    EXPECT_NEAR(LegendrePolynomial(4, x), expected, 1e-15);
  }
}

TEST(LegendrePolynomial, P5ClosedForm) {
  for (double x : kAllX) {
    double x2 = x * x;
    double x3 = x2 * x;
    double x5 = x3 * x2;
    double expected = (63.0 * x5 - 70.0 * x3 + 15.0 * x) / 8.0;
    EXPECT_NEAR(LegendrePolynomial(5, x), expected, 1e-15);
  }
}

TEST(LegendrePolynomial, EndpointValues) {
  for (int n = 0; n <= 10; ++n) {
    EXPECT_DOUBLE_EQ(LegendrePolynomial(n, 1.0), 1.0) << "n=" << n;
    double expected_minus_one = (n % 2 == 0) ? 1.0 : -1.0;
    EXPECT_DOUBLE_EQ(LegendrePolynomial(n, -1.0), expected_minus_one)
        << "n=" << n;
  }
}

// Bonnet's recursion, as an independent sanity check of the recursive
// implementation across orders: (n+1) P_{n+1} = (2n+1) x P_n - n P_{n-1}
TEST(LegendrePolynomial, BonnetRecurrence) {
  for (double x : kAllX) {
    for (int n = 1; n <= 9; ++n) {
      double lhs = (n + 1) * LegendrePolynomial(n + 1, x);
      double rhs = (2 * n + 1) * x * LegendrePolynomial(n, x) -
                   n * LegendrePolynomial(n - 1, x);
      EXPECT_DOUBLE_EQ(lhs, rhs);
    }
  }
}

TEST(LegendrePolynomial, ParitySymmetry) {
  for (double x : kAllX) {
    for (int n = 0; n <= 8; ++n) {
      double sign = (n % 2 == 0) ? 1.0 : -1.0;
      EXPECT_DOUBLE_EQ(LegendrePolynomial(n, -x),
                       sign * LegendrePolynomial(n, x));
    }
  }
}

// =======================================================================
// LegendrePolynomialDerivative(n, x)
// =======================================================================

TEST(LegendrePolynomialDerivative, D0IsZero) {
  for (double x : kAllX) {
    EXPECT_DOUBLE_EQ(LegendrePolynomialDerivative(0, x), 0.0);
  }
}

// P_1(x) = x, so P_1'(x) should be identically 1 everywhere, including at
// x = 0.5 where 1/(1-x^2) = 4/3 != 1. This test is expected to fail
// against the current implementation, which hardcodes the n==1 case as
// 1.0 / (1.0 - x * x) instead of falling through to the general formula
// (which itself simplifies to 1 for n=1).
TEST(LegendrePolynomialDerivative, D1IsOne) {
  for (double x : kInteriorX) {
    EXPECT_DOUBLE_EQ(LegendrePolynomialDerivative(1, x), 1.0);
  }
}

TEST(LegendrePolynomialDerivative, D2ClosedForm) {
  for (double x : kInteriorX) {
    EXPECT_DOUBLE_EQ(LegendrePolynomialDerivative(2, x), 3.0 * x);
  }
}

TEST(LegendrePolynomialDerivative, D3ClosedForm) {
  for (double x : kInteriorX) {
    double expected = 0.5 * (15.0 * x * x - 3.0);
    EXPECT_DOUBLE_EQ(LegendrePolynomialDerivative(3, x), expected);
  }
}

TEST(LegendrePolynomialDerivative, D4ClosedForm) {
  for (double x : kInteriorX) {
    double expected = 0.5 * (35.0 * x * x * x - 15.0 * x);
    EXPECT_DOUBLE_EQ(LegendrePolynomialDerivative(4, x), expected);
  }
}

// Cross-check the derivative against central finite differences of
// LegendrePolynomial itself. Restricted to the interior since the
// closed-form recurrence divides by (1 - x^2), which is singular at the
// endpoints.
TEST(LegendrePolynomialDerivative, MatchesFiniteDifference) {
  const double h = 1e-6;
  for (int n = 0; n <= 8; ++n) {
    for (double x : kInteriorX) {
      double numerical =
          (LegendrePolynomial(n, x + h) - LegendrePolynomial(n, x - h)) /
          (2.0 * h);
      double analytic = LegendrePolynomialDerivative(n, x);
      EXPECT_NEAR(analytic, numerical, 1e-9);
    }
  }
}

// Legendre's differential equation:
//   (1 - x^2) P_n''(x) - 2x P_n'(x) + n(n+1) P_n(x) = 0
// P_n'' is approximated via a central finite difference of
// LegendrePolynomialDerivative, so this also cross-validates the two
// functions against each other through an independent identity.
TEST(LegendrePolynomialDerivative, SatisfiesLegendreODE) {
  const double h = 1e-4;
  for (int n = 0; n <= 6; ++n) {
    for (double x : kInteriorX) {
      double d1 = LegendrePolynomialDerivative(n, x);
      double d2 =
          (LegendrePolynomial(n, x + h) - 2.0 * LegendrePolynomial(n, x) +
           LegendrePolynomial(n, x - h)) /
          (h * h);
      double lhs = (1.0 - x * x) * d2 - 2.0 * x * d1 +
                   n * (n + 1) * LegendrePolynomial(n, x);
      EXPECT_NEAR(lhs, 0.0, 1e-6);
    }
  }
}

// Endpoint derivative values: P_n'(1) = n(n+1)/2, P_n'(-1) = (-1)^(n+1) *
// n(n+1)/2.
TEST(LegendrePolynomialDerivative, EndpointDerivativeAtPlusOne) {
  for (int n = 0; n <= 10; ++n) {
    double expected = n * (n + 1) / 2.0;
    EXPECT_DOUBLE_EQ(LegendrePolynomialDerivative(n, 1.0), expected)
        << "n=" << n;
  }
}

TEST(LegendrePolynomialDerivative, EndpointDerivativeAtMinusOne) {
  for (int n = 0; n <= 10; ++n) {
    double sign = (n % 2 == 0) ? -1.0 : 1.0;
    double expected = sign * n * (n + 1) / 2.0;
    EXPECT_DOUBLE_EQ(LegendrePolynomialDerivative(n, -1.0), expected)
        << "n=" << n;
  }
}

}  // namespace hummingbird::math
