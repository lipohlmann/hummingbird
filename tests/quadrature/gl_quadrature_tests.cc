#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <vector>

#include "quadrature/gl_quadrature.h"
#include "quadrature/quadrature.h"

namespace hummingbird::quadrature {

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Builds the QuadraturePair vector required by Integrate() by evaluating
// `f` at every abscissa of `quad`.
std::vector<QuadraturePair> EvaluateAt(const GLQuadrature& quad,
                                       const std::function<double(double)>& f) {
  std::vector<QuadraturePair> pairs;
  pairs.reserve(quad.abscissas().size());
  for (auto i = 0; i < quad.abscissas().size(); i++) {
    pairs.push_back(QuadraturePair{i, f(quad.GetAbscissa(i))});
  }
  return pairs;
}

// Closed-form value of \int_{-1}^{1} x^k dx.
double AnalyticMonomialIntegral(const unsigned int k) {
  if (k % 2 == 1) {
    return 0.0;
  }
  return 2.0 / static_cast<double>(k + 1);
}

// ---------------------------------------------------------------------------
// Parameterized tests: exactness on monomials / polynomials
// ---------------------------------------------------------------------------

// Parameter is the number of quadrature points, n.
class GLQuadratureExactnessTest
    : public ::testing::TestWithParam<unsigned int> {};

TEST_P(GLQuadratureExactnessTest, IntegratesAllMonomialsUpToDegree2NMinus1) {
  const unsigned int n = GetParam();
  GLQuadrature quad(n);

  ASSERT_EQ(quad.abscissas().size(), n);

  // An n-point Gauss-Legendre rule is exact for all polynomials of degree
  // <= 2n - 1, so every monomial x^k in that range must integrate exactly.
  const unsigned int max_exact_degree = 2 * n - 1;
  for (unsigned int k = 0; k <= max_exact_degree; ++k) {
    auto pairs = EvaluateAt(quad, [k](double x) { return std::pow(x, k); });
    const double result = quad.Integrate(pairs);
    const double expected = AnalyticMonomialIntegral(k);
    EXPECT_DOUBLE_EQ(result, expected)
        << "Failed for n=" << n << ", monomial degree k=" << k;
  }
}

TEST_P(GLQuadratureExactnessTest, IntegratesArbitraryPolynomialInRange) {
  const unsigned int n = GetParam();
  GLQuadrature quad(n);

  // Build a polynomial that uses every degree from 0 up to the max exact
  // degree, e.g. p(x) = sum_{k=0}^{2n-1} (k + 1) * x^k
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
    expected += static_cast<double>(k + 1) * AnalyticMonomialIntegral(k);
  }

  auto pairs = EvaluateAt(quad, poly);
  const double result = quad.Integrate(pairs);
  EXPECT_DOUBLE_EQ(result, expected);
}

TEST_P(GLQuadratureExactnessTest, WeightsSumToTwo) {
  // \int_{-1}^{1} 1 dx = 2, i.e. sum of all weights must equal 2.
  const unsigned int n = GetParam();
  GLQuadrature quad(n);

  auto pairs = EvaluateAt(quad, [](double) { return 1.0; });
  const double result = quad.Integrate(pairs);
  EXPECT_DOUBLE_EQ(result, 2.0);
}

TEST_P(GLQuadratureExactnessTest, WeightsArePositive) {
  const unsigned int n = GetParam();
  GLQuadrature quad(n);

  for (const double x : quad.abscissas()) {
    EXPECT_GT(quad.GetWeight(x), 0.0)
        << "Gauss-Legendre weights must all be strictly positive.";
  }
}

TEST_P(GLQuadratureExactnessTest, AbscissasLieWithinIntervalAndAreSymmetric) {
  const unsigned int n = GetParam();
  GLQuadrature quad(n);
  const auto& abscissas = quad.abscissas();

  std::vector<double> sorted(abscissas.begin(), abscissas.end());
  std::sort(sorted.begin(), sorted.end());

  for (double x : sorted) {
    EXPECT_GE(x, -1.0);
    EXPECT_LE(x, 1.0);
  }

  // Gauss-Legendre nodes are symmetric about 0: x_i == -x_{n-1-i}.
  for (std::size_t i = 0; i < sorted.size(); ++i) {
    const double mirrored = sorted[sorted.size() - 1 - i];
    EXPECT_NEAR(sorted[i], -mirrored, 1e-12);
  }
}

INSTANTIATE_TEST_SUITE_P(VariousPointCounts, GLQuadratureExactnessTest,
                         ::testing::Values(1u, 2u, 3u, 4u, 5u, 6u, 8u, 10u));

// ---------------------------------------------------------------------------
// Boundary / negative test: a polynomial one degree too high should NOT
// generally be integrated exactly.
// ---------------------------------------------------------------------------

TEST(GLQuadratureBoundaryTest,
     DoesNotExactlyIntegratePolynomialAboveDegreeBound) {
  const unsigned int n = 3;
  GLQuadrature quad(n);

  // Degree 2n = 6 is one past the exactness bound (2n - 1 = 5).
  const unsigned int degree = 2 * n;
  auto pairs =
      EvaluateAt(quad, [degree](double x) { return std::pow(x, degree); });
  const double result = quad.Integrate(pairs);
  const double expected = AnalyticMonomialIntegral(degree);

  // We expect a real discrepancy here, not just floating point noise.
  EXPECT_GT(std::abs(result - expected), 1e-3);
}

// ---------------------------------------------------------------------------
// Non-polynomial test: exp(x) has a known closed form and provides a
// meaningful accuracy/convergence check for a transcendental integrand.
// ---------------------------------------------------------------------------

TEST(GLQuadratureNonPolynomialTest,
     ApproximatesExponentialIntegralWithinTolerance) {
  // \int_{-1}^{1} e^x dx = e - 1/e
  const double expected = std::exp(1.0) - std::exp(-1.0);

  const unsigned int n = 20;  // enough points for e^x to converge tightly
  GLQuadrature quad(n);

  auto pairs = EvaluateAt(quad, [](double x) { return std::exp(x); });
  const double result = quad.Integrate(pairs);

  EXPECT_NEAR(result, expected, 1e-10);
}

TEST(GLQuadratureNonPolynomialTest,
     ExponentialIntegralConvergesAsPointsIncrease) {
  const double expected = std::exp(1.0) - std::exp(-1.0);

  double previous_error = std::numeric_limits<double>::max();
  for (unsigned int n : {2u, 4u, 6u, 8u, 10u}) {
    GLQuadrature quad(n);
    auto pairs = EvaluateAt(quad, [](double x) { return std::exp(x); });
    const double result = quad.Integrate(pairs);
    const double error = std::abs(result - expected);

    EXPECT_LE(error, previous_error)
        << "Error should not increase as n grows (n=" << n << ").";
    previous_error = error;
  }

  // By n=10 the error should be essentially machine precision.
  EXPECT_LT(previous_error, 1e-12);
}

}  // namespace hummingbird::quadrature