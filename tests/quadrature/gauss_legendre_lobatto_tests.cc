#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

#include "../tests/quadrature/helpers.h"
#include "quadrature/gauss_legendre_lobatto.h"
#include "quadrature/quadrature_base.h"
#include "utils/constants.h"

namespace hummingbird::quadrature {

// ---------------------------------------------------------------------------
// Parameterized tests: exactness on monomials / polynomials
//
// An N-point Gauss-Legendre-Lobatto rule fixes two abscissas at the
// endpoints, so it is exact only up to degree 2N - 3 (two degrees of
// freedom fewer than the (2N - 1)-exact Gauss-Legendre rule).
// ---------------------------------------------------------------------------

// Parameter is the number of quadrature points, n. n must be >= 2 since GLL
// requires at least the two endpoints.
class GLLQuadratureExactnessTest
    : public ::testing::TestWithParam<unsigned int> {};

TEST_P(GLLQuadratureExactnessTest, IntegratesAllMonomialsUpToDegree2NMinus3) {
  const unsigned int n = GetParam();
  GaussLegendreLobatto quad(n);

  ASSERT_EQ(quad.abscissas().size(), n);

  const int max_exact_degree = 2 * static_cast<int>(n) - 3;
  for (int k = 0; k <= max_exact_degree; ++k) {
    auto pairs = GLLEvaluateAt(quad, [k](double x) { return std::pow(x, k); });
    const double result = quad.Integrate(pairs);
    const double expected = AnalyticMonomialIntegral(k);
    EXPECT_NEAR(result, expected, utils::TOLERANCE)
        << "Failed for n=" << n << ", monomial degree k=" << k;
  }
}

TEST_P(GLLQuadratureExactnessTest, IntegratesArbitraryPolynomialInRange) {
  const unsigned int n = GetParam();
  GaussLegendreLobatto quad(n);

  const int max_exact_degree = 2 * static_cast<int>(n) - 3;
  auto poly = [max_exact_degree](double x) {
    double sum = 0.0;
    for (int k = 0; k <= max_exact_degree; ++k) {
      sum += static_cast<double>(k + 1) * std::pow(x, k);
    }
    return sum;
  };

  double expected = 0.0;
  for (int k = 0; k <= max_exact_degree; ++k) {
    expected += static_cast<double>(k + 1) * AnalyticMonomialIntegral(k);
  }

  auto pairs = GLLEvaluateAt(quad, poly);
  const double result = quad.Integrate(pairs);
  EXPECT_DOUBLE_EQ(result, expected);
}

TEST_P(GLLQuadratureExactnessTest, WeightsSumToTwo) {
  // \int_{-1}^{1} 1 dx = 2, i.e. sum of all weights must equal 2.
  const unsigned int n = GetParam();
  GaussLegendreLobatto quad(n);

  auto pairs = GLLEvaluateAt(quad, [](double) { return 1.0; });
  const double result = quad.Integrate(pairs);
  EXPECT_DOUBLE_EQ(result, 2.0);
}

TEST_P(GLLQuadratureExactnessTest, WeightsArePositive) {
  const unsigned int n = GetParam();
  GaussLegendreLobatto quad(n);

  for (unsigned int i = 0; i < quad.abscissas().size(); ++i) {
    EXPECT_GT(quad.GetWeight(i), 0.0)
        << "Gauss-Legendre-Lobatto weights must all be strictly positive "
           "(index "
        << i << ").";
  }
}

TEST_P(GLLQuadratureExactnessTest, AbscissasLieWithinIntervalAndAreSymmetric) {
  const unsigned int n = GetParam();
  GaussLegendreLobatto quad(n);
  const auto& abscissas = quad.abscissas();

  std::vector<double> sorted(abscissas.begin(), abscissas.end());
  std::sort(sorted.begin(), sorted.end());

  for (double x : sorted) {
    EXPECT_GE(x, -1.0);
    EXPECT_LE(x, 1.0);
  }

  // GLL nodes are symmetric about 0: x_i == -x_{n-1-i}.
  for (std::size_t i = 0; i < sorted.size(); ++i) {
    const double mirrored = sorted[sorted.size() - 1 - i];
    EXPECT_DOUBLE_EQ(sorted[i], -mirrored);
  }
}

TEST_P(GLLQuadratureExactnessTest, EndpointsAreFixedAtPlusMinusOne) {
  // Unlike Gauss-Legendre, GLL always includes both endpoints of [-1, 1]
  // as abscissas.
  const unsigned int n = GetParam();
  GaussLegendreLobatto quad(n);
  const auto& abscissas = quad.abscissas();

  const double min_val = abscissas.front();
  const double max_val = abscissas.back();

  EXPECT_DOUBLE_EQ(min_val, -1.0);
  EXPECT_DOUBLE_EQ(max_val, 1.0);
}

INSTANTIATE_TEST_SUITE_P(VariousPointCounts, GLLQuadratureExactnessTest,
                         ::testing::Values(2u, 3u, 4u, 5u, 6u, 8u, 10u));

// ---------------------------------------------------------------------------
// Two-point edge case: with only the two endpoints, GLL degenerates into
// the trapezoidal rule (weights of 1 each), exact only for linear functions.
// ---------------------------------------------------------------------------

TEST(GLLQuadratureTwoPointTest, DegeneratesToTrapezoidalRule) {
  GaussLegendreLobatto quad(2);
  ASSERT_EQ(quad.abscissas().size(), 2u);

  for (unsigned int i = 0; i < 2; ++i) {
    EXPECT_DOUBLE_EQ(quad.GetWeight(i), 1.0);
  }
}

// ---------------------------------------------------------------------------
// Boundary / negative test: a polynomial one degree too high should NOT
// generally be integrated exactly.
// ---------------------------------------------------------------------------

TEST(GLLQuadratureBoundaryTest,
     DoesNotExactlyIntegratePolynomialAboveDegreeBound) {
  const unsigned int n = 4;
  GaussLegendreLobatto quad(n);

  // Exactness bound is 2n - 3 = 5; degree 2n - 2 = 6 is one past it.
  const unsigned int degree = 2 * n - 2;
  auto pairs =
      GLLEvaluateAt(quad, [degree](double x) { return std::pow(x, degree); });
  const double result = quad.Integrate(pairs);
  const double expected = AnalyticMonomialIntegral(degree);

  // We expect a real discrepancy here, not just floating point noise.
  EXPECT_GT(std::abs(result - expected), 1e-3);
}

// ---------------------------------------------------------------------------
// Non-polynomial test: exp(x) has a known closed form and provides a
// meaningful accuracy/convergence check for a transcendental integrand.
// ---------------------------------------------------------------------------

TEST(GLLQuadratureNonPolynomialTest,
     ApproximatesExponentialIntegralWithinTolerance) {
  // \int_{-1}^{1} e^x dx = e - 1/e
  const double expected = std::exp(1.0) - std::exp(-1.0);

  const unsigned int n = 20;  // enough points for e^x to converge tightly
  GaussLegendreLobatto quad(n);

  auto pairs = GLLEvaluateAt(quad, [](double x) { return std::exp(x); });
  const double result = quad.Integrate(pairs);

  EXPECT_DOUBLE_EQ(result, expected);
}

TEST(GLLQuadratureNonPolynomialTest,
     ExponentialIntegralConvergesAsPointsIncrease) {
  const double expected = std::exp(1.0) - std::exp(-1.0);

  double previous_error = std::numeric_limits<double>::max();
  for (unsigned int n : {3u, 5u, 7u, 9u}) {
    GaussLegendreLobatto quad(n);
    auto pairs = GLLEvaluateAt(quad, [](double x) { return std::exp(x); });
    const double result = quad.Integrate(pairs);
    const double error = std::abs(result - expected);

    EXPECT_LE(error, previous_error)
        << "Error should not increase as n grows (n=" << n << ").";
    previous_error = error;
  }

  // By n=9 the error should be very small (GLL converges slightly slower
  // than Gauss-Legendre for the same n, since two degrees of freedom are
  // spent fixing the endpoints).
  EXPECT_LT(previous_error, 1e-8);
}

}  // namespace hummingbird::quadrature
