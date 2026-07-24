#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

#include "../tests/quadrature/helpers.h"
#include "quadrature/gauss_chebyshev.h"
#include "quadrature/quadrature_base.h"
#include "utils/constants.h"

namespace hummingbird::quadrature {

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Product form of the central binomial coefficient C(2m, m), computed
// iteratively to avoid overflow: C(2m, m) = prod_{i=1}^{m} (m + i) / i.
double CentralBinomialCoefficient(const unsigned int m) {
  double result = 1.0;
  for (unsigned int i = 1; i <= m; ++i) {
    result *= static_cast<double>(m + i) / static_cast<double>(i);
  }
  return result;
}

// Closed-form value of \int_{-1}^{1} x^k / sqrt(1 - x^2) dx.
// Odd powers vanish by symmetry. For even k = 2m:
//   \int_{-1}^{1} x^{2m} / sqrt(1 - x^2) dx = pi * C(2m, m) / 4^m
double AnalyticWeightedMonomialIntegral(const unsigned int k) {
  if (k % 2 == 1) {
    return 0.0;
  }
  const unsigned int m = k / 2;
  return M_PI * CentralBinomialCoefficient(m) / std::pow(4.0, m);
}

// ---------------------------------------------------------------------------
// Parameterized tests: exactness on monomials / polynomials
//
// Gauss-Chebyshev is a genuine N-point Gaussian rule, so like
// Gauss-Legendre it is exact for polynomials up to degree 2N - 1, measured
// against the weighted integral above.
// ---------------------------------------------------------------------------

// Parameter is the number of quadrature points, n.
class GCQuadratureExactnessTest
    : public ::testing::TestWithParam<unsigned int> {};

TEST_P(GCQuadratureExactnessTest, IntegratesAllMonomialsUpToDegree2NMinus1) {
  const unsigned int n = GetParam();
  GaussChebyshev quad(n);

  ASSERT_EQ(quad.abscissas().size(), n);

  const unsigned int max_exact_degree = 2 * n - 1;
  for (unsigned int k = 0; k <= max_exact_degree; ++k) {
    auto pairs = EvaluateAt(quad, [k](double x) { return std::pow(x, k); });
    const double result = quad.Integrate(pairs);
    const double expected = AnalyticWeightedMonomialIntegral(k);
    EXPECT_NEAR(result, expected, utils::TOLERANCE)
        << "Failed for n=" << n << ", monomial degree k=" << k;
  }
}

TEST_P(GCQuadratureExactnessTest, IntegratesArbitraryPolynomialInRange) {
  const unsigned int n = GetParam();
  GaussChebyshev quad(n);

  const unsigned int max_exact_degree = 2 * n - 1;
  auto poly = [max_exact_degree](double x) {
    double sum = 0.0;
    for (unsigned int k = 0; k <= max_exact_degree; ++k) {
      sum += static_cast<double>(k + 1) * std::pow(x, k);
    }
    return sum;
  };

  double expected = 0.0;
  for (unsigned int k = 0; k <= max_exact_degree; ++k) {
    expected +=
        static_cast<double>(k + 1) * AnalyticWeightedMonomialIntegral(k);
  }

  auto pairs = EvaluateAt(quad, poly);
  const double result = quad.Integrate(pairs);
  EXPECT_DOUBLE_EQ(result, expected);
}

TEST_P(GCQuadratureExactnessTest, WeightsSumToPi) {
  // \int_{-1}^{1} 1 / sqrt(1 - x^2) dx = pi.
  const unsigned int n = GetParam();
  GaussChebyshev quad(n);

  auto pairs = EvaluateAt(quad, [](double) { return 1.0; });
  const double result = quad.Integrate(pairs);
  EXPECT_DOUBLE_EQ(result, M_PI);
}

TEST_P(GCQuadratureExactnessTest, AllWeightsAreEqualToPiOverN) {
  // Unlike GL/GLL, every Gauss-Chebyshev weight is identical: w_k = pi / n.
  const unsigned int n = GetParam();
  GaussChebyshev quad(n);

  const double expected_weight = M_PI / static_cast<double>(n);
  for (unsigned int i = 0; i < n; ++i) {
    EXPECT_DOUBLE_EQ(quad.GetWeight(i), expected_weight)
        << "Weight mismatch at index " << i << " for n=" << n;
  }
}

TEST_P(GCQuadratureExactnessTest, WeightsArePositive) {
  const unsigned int n = GetParam();
  GaussChebyshev quad(n);

  for (unsigned int i = 0; i < n; ++i) {
    EXPECT_GT(quad.GetWeight(i), 0.0)
        << "Gauss-Chebyshev weights must all be strictly positive (index " << i
        << ").";
  }
}

TEST_P(GCQuadratureExactnessTest,
       AbscissasLieStrictlyWithinIntervalAndAreSymmetric) {
  const unsigned int n = GetParam();
  GaussChebyshev quad(n);
  const auto& abscissas = quad.abscissas();

  std::vector<double> sorted(abscissas.begin(), abscissas.end());
  std::sort(sorted.begin(), sorted.end());

  for (double x : sorted) {
    // Chebyshev nodes cos((2k+1)pi/2N) never reach the endpoints exactly,
    // unlike Gauss-Legendre-Lobatto.
    EXPECT_GT(x, -1.0);
    EXPECT_LT(x, 1.0);
  }

  // Gauss-Chebyshev nodes are symmetric about 0: x_i == -x_{n-1-i}.
  for (std::size_t i = 0; i < sorted.size(); ++i) {
    const double mirrored = sorted[sorted.size() - 1 - i];
    EXPECT_NEAR(sorted[i], -mirrored, utils::TOLERANCE);
  }
}

INSTANTIATE_TEST_SUITE_P(VariousPointCounts, GCQuadratureExactnessTest,
                         ::testing::Values(1u, 2u, 3u, 4u, 5u, 6u, 8u, 10u));

// ---------------------------------------------------------------------------
// Known closed-form abscissa check for a small, hand-verifiable case.
// ---------------------------------------------------------------------------

TEST(GCQuadratureExplicitFormulaTest, MatchesKnownThreePointAbscissas) {
  // xi_k = cos((2k+1)pi / (2N)) for N = 3:
  //   xi_0 = cos(pi/6)  =  sqrt(3)/2
  //   xi_1 = cos(pi/2)  =  0
  //   xi_2 = cos(5pi/6) = -sqrt(3)/2
  const unsigned int n = 3;
  GaussChebyshev quad(n);

  std::vector<double> sorted(quad.abscissas().begin(), quad.abscissas().end());
  std::sort(sorted.begin(), sorted.end());

  EXPECT_DOUBLE_EQ(sorted[0], -std::sqrt(3.0) / 2.0);
  EXPECT_NEAR(sorted[1], 0.0, utils::TOLERANCE);
  EXPECT_DOUBLE_EQ(sorted[2], std::sqrt(3.0) / 2.0);
}

// ---------------------------------------------------------------------------
// Boundary / negative test: a polynomial one degree too high should NOT
// generally be integrated exactly.
// ---------------------------------------------------------------------------

TEST(GCQuadratureBoundaryTest,
     DoesNotExactlyIntegratePolynomialAboveDegreeBound) {
  const unsigned int n = 3;
  GaussChebyshev quad(n);

  // Degree 2n = 6 is one past the exactness bound (2n - 1 = 5).
  const unsigned int degree = 2 * n;
  auto pairs =
      EvaluateAt(quad, [degree](double x) { return std::pow(x, degree); });
  const double result = quad.Integrate(pairs);
  const double expected = AnalyticWeightedMonomialIntegral(degree);

  // We expect a real discrepancy here, not just floating point noise.
  EXPECT_GT(std::abs(result - expected), 1e-3);
}

// ---------------------------------------------------------------------------
// Non-polynomial test: exp(x) has a known closed form (a modified Bessel
// function) and provides a meaningful accuracy/convergence check for a
// transcendental integrand weighted by 1/sqrt(1-x^2).
// ---------------------------------------------------------------------------

TEST(GCQuadratureNonPolynomialTest,
     ExponentialIntegralConvergesAsPointsIncrease) {
  // \int_{-1}^{1} e^x / sqrt(1 - x^2) dx = pi * I_0(1), where I_0 is the
  // modified Bessel function of the first kind, order 0.
  // I_0(1) ~= 1.2660658777520084
  const double expected = M_PI * 1.2660658777520084;

  double previous_error = std::numeric_limits<double>::max();
  for (unsigned int n : {2u, 4u, 6u, 8u, 10u}) {
    GaussChebyshev quad(n);
    auto pairs = EvaluateAt(quad, [](double x) { return std::exp(x); });
    const double result = quad.Integrate(pairs);
    const double error = std::abs(result - expected);

    EXPECT_LE(error, previous_error)
        << "Error should not increase as n grows (n=" << n << ").";
    previous_error = error;
  }

  // By n=10 the error should be essentially machine precision.
  EXPECT_LT(previous_error, 1e-10);
}

}  // namespace hummingbird::quadrature