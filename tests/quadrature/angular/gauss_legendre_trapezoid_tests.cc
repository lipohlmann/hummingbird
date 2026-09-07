#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <limits>
#include <tuple>
#include <vector>

#include "quadrature/angular/gauss_legendre_trapezoid.h"
#include "quadrature/angular/ordinate.h"
#include "quadrature/gauss_legendre.h"
#include "utils/constants.h"

namespace hummingbird {

namespace {

constexpr double kTwoPi = 2.0 * M_PI;
constexpr double kFourPi = 4.0 * M_PI;

// Sums a per-ordinate quantity over every point in the quadrature, weighted
// by that point's quadrature weight. This plays the same role as
// QuadratureBase::Integrate(), but works directly off direction cosines
// since Ordinate (unlike double) isn't a scalar abscissa type the shared
// helper machinery expects.
double WeightedSum(const GaussLegendreTrapezoid& quad,
                   const std::function<double(Ordinate)>& f) {
  double sum = 0.0;
  for (unsigned int i = 0; i < quad.n_points(); ++i) {
    Ordinate ord = quad.GetAbscissa(i);
    sum += quad.GetWeight(i) * f(ord);
  }
  return sum;
}

// Matches the implementation's flattening of (polar index i, azimuthal
// index j) into a single point index: the polar loop is outer, so each
// polar level owns a contiguous run of n_azim azimuthal points.
size_t Index(size_t n_azim, size_t i, size_t j) { return i * n_azim + j; }

// Azimuthal angle in [0, 2*pi), recovered from the ordinate itself (as
// opposed to the construction formula) so tests using this are genuinely
// checking the quadrature's output.
double Phi(Ordinate ord) {
  double phi = std::atan2(ord.y(), ord.x());
  if (phi < 0.0) phi += kTwoPi;
  return phi;
}

// n!! with the convention (-1)!! = 1, used by the closed-form sphere-monomial
// integral below.
double DoubleFactorial(int n) {
  double result = 1.0;
  for (int i = n; i > 0; i -= 2) result *= i;
  return result;
}

// Closed-form value of int_{S^2} x^a y^b z^c dOmega. Vanishes unless a, b, c
// are all even, in which case it equals
// 4*pi * (a-1)!! * (b-1)!! * (c-1)!! / (a+b+c+1)!!
// This is a standard result for moments of monomials over the unit sphere
// and gives an independent, closed-form target to check the quadrature
// against (as opposed to re-deriving the expected value from the
// quadrature's own weights).
double AnalyticSphereMonomialIntegral(int a, int b, int c) {
  if (a % 2 != 0 || b % 2 != 0 || c % 2 != 0) return 0.0;
  return 4.0 * M_PI * DoubleFactorial(a - 1) * DoubleFactorial(b - 1) *
         DoubleFactorial(c - 1) / DoubleFactorial(a + b + c + 1);
}

struct Monomial {
  int a, b, c;
};

}  // namespace

// ---------------------------------------------------------------------------
// Structural tests
// ---------------------------------------------------------------------------

class GLTStructureTest
    : public ::testing::TestWithParam<std::tuple<size_t, size_t>> {};

TEST_P(GLTStructureTest, ProducesExpectedNumberOfOrdinates) {
  auto [n_azim, n_polar] = GetParam();
  GaussLegendreTrapezoid quad(n_azim, n_polar);
  EXPECT_EQ(quad.n_points(), n_azim * n_polar);
}

TEST_P(GLTStructureTest, EveryOrdinateIsAUnitVector) {
  auto [n_azim, n_polar] = GetParam();
  GaussLegendreTrapezoid quad(n_azim, n_polar);
  for (unsigned int i = 0; i < quad.n_points(); ++i) {
    Ordinate ord = quad.GetAbscissa(i);
    const double norm_sq =
        ord.x() * ord.x() + ord.y() * ord.y() + ord.z() * ord.z();
    EXPECT_NEAR(norm_sq, 1.0, EXP_NEAR_TOLERANCE)
        << "Direction cosines at index " << i << " are not normalized.";
  }
}

TEST_P(GLTStructureTest, AllWeightsArePositive) {
  auto [n_azim, n_polar] = GetParam();
  GaussLegendreTrapezoid quad(n_azim, n_polar);
  for (unsigned int i = 0; i < quad.n_points(); ++i) {
    EXPECT_GT(quad.GetWeight(i), 0.0) << "Non-positive weight at index " << i;
  }
}

TEST_P(GLTStructureTest, WeightsSumToFourPi) {
  // The quadrature approximates the surface integral over the unit sphere,
  // whose total solid angle is 4*pi. This holds for any n_azim/n_polar
  // combination, since it only relies on the Gauss-Legendre weights summing
  // to 2 (length of [-1,1]) and the trapezoid weights summing to 2*pi (full
  // azimuthal circle).
  auto [n_azim, n_polar] = GetParam();
  GaussLegendreTrapezoid quad(n_azim, n_polar);
  const double sum = WeightedSum(quad, [](Ordinate) { return 1.0; });
  EXPECT_NEAR(sum, kFourPi, EXP_NEAR_TOLERANCE * 10.0);
}

TEST_P(GLTStructureTest, IntegrateAgreesWithManualWeightedSumForAConstant) {
  // Exercise the inherited QuadratureBase::Integrate() / QuadraturePair path
  // directly, since that's the interface production code will actually use.
  auto [n_azim, n_polar] = GetParam();
  GaussLegendreTrapezoid quad(n_azim, n_polar);

  std::vector<QuadraturePair> pairs;
  pairs.reserve(quad.n_points());
  for (unsigned int i = 0; i < quad.n_points(); ++i) {
    pairs.push_back({static_cast<int>(i), 1.0});
  }

  EXPECT_NEAR(quad.Integrate(pairs), kFourPi, EXP_NEAR_TOLERANCE);
}

TEST_P(GLTStructureTest,
       WeightsMatchProductOfPolarWeightAndUniformAzimuthalWeight) {
  // weight(i, j) == w_polar(i) * (2*pi / n_azim) for every azimuthal index j
  // at a given polar level i, straight from the constructor.
  auto [n_azim, n_polar] = GetParam();
  GaussLegendreTrapezoid quad(n_azim, n_polar);
  GaussLegendre gl_polar(n_polar);
  const double delta_azim = kTwoPi / n_azim;

  for (size_t i = 0; i < n_polar; ++i) {
    const double expected_weight = gl_polar.GetWeight(i) * delta_azim;
    for (size_t j = 0; j < n_azim; ++j) {
      EXPECT_NEAR(quad.GetWeight(Index(n_azim, i, j)), expected_weight,
                  EXP_NEAR_TOLERANCE)
          << "Mismatch at polar index " << i << ", azimuthal index " << j;
    }
  }
}

TEST_P(GLTStructureTest, OrdinatesMatchClosedFormPolarAndAzimuthalAngles) {
  // Each ordinate should be (sin(theta_i)*cos(phi_j), sin(theta_i)*sin(phi_j),
  // mu_i), with mu_i the i-th Gauss-Legendre node, theta_i = acos(mu_i), and
  // phi_j = -pi + j * (2*pi / n_azim).
  auto [n_azim, n_polar] = GetParam();
  GaussLegendreTrapezoid quad(n_azim, n_polar);
  GaussLegendre gl_polar(n_polar);
  const double delta_azim = kTwoPi / n_azim;

  for (size_t i = 0; i < n_polar; ++i) {
    const double mu = gl_polar.GetAbscissa(i);
    const double theta = std::acos(mu);
    for (size_t j = 0; j < n_azim; ++j) {
      const double phi = -M_PI + static_cast<double>(j) * delta_azim;
      Ordinate ord = quad.GetAbscissa(Index(n_azim, i, j));
      EXPECT_NEAR(ord.x(), std::sin(theta) * std::cos(phi), EXP_NEAR_TOLERANCE)
          << "at polar index " << i << ", azimuthal index " << j;
      EXPECT_NEAR(ord.y(), std::sin(theta) * std::sin(phi), EXP_NEAR_TOLERANCE)
          << "at polar index " << i << ", azimuthal index " << j;
      EXPECT_NEAR(ord.z(), mu, EXP_NEAR_TOLERANCE)
          << "at polar index " << i << ", azimuthal index " << j;
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    VariousOrders, GLTStructureTest,
    ::testing::Values(std::make_tuple(1, 1), std::make_tuple(1, 2),
                      std::make_tuple(2, 1), std::make_tuple(2, 2),
                      std::make_tuple(3, 4), std::make_tuple(4, 3),
                      std::make_tuple(5, 5)));

// ---------------------------------------------------------------------------
// Moment tests
//
// For a well-formed quadrature over the sphere, the first moments of the
// direction cosines and all cross moments must vanish by symmetry, and the
// second "diagonal" moments int(Omega_x^2), int(Omega_y^2), int(Omega_z^2)
// must each equal 4*pi/3 -- provided there are enough points to avoid
// aliasing. The periodic trapezoid rule with n_azim equally spaced points
// integrates exp(i*k*phi) exactly unless n_azim divides k, so first moments
// (k=1) need n_azim >= 2, and the diagonal second moments (which carry a
// cos(2*phi)/sin(2*phi) harmonic, k=2) need n_azim >= 3. n_polar >= 2 is
// needed for the polar moment, since mu^2 is a degree-2 polynomial.
// ---------------------------------------------------------------------------

class GLTMomentTest
    : public ::testing::TestWithParam<std::tuple<size_t, size_t>> {};

TEST_P(GLTMomentTest, FirstMomentsVanish) {
  auto [n_azim, n_polar] = GetParam();
  GaussLegendreTrapezoid quad(n_azim, n_polar);

  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.x(); }), 0.0,
              EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.y(); }), 0.0,
              EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.z(); }), 0.0,
              EXP_NEAR_TOLERANCE);
}

TEST_P(GLTMomentTest, CrossMomentsVanish) {
  auto [n_azim, n_polar] = GetParam();
  GaussLegendreTrapezoid quad(n_azim, n_polar);

  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.x() * o.y(); }), 0.0,
              EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.y() * o.z(); }), 0.0,
              EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.x() * o.z(); }), 0.0,
              EXP_NEAR_TOLERANCE);
}

TEST_P(GLTMomentTest, DiagonalSecondMomentsEqualFourPiOverThree) {
  auto [n_azim, n_polar] = GetParam();
  GaussLegendreTrapezoid quad(n_azim, n_polar);
  const double expected = kFourPi / 3.0;

  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.x() * o.x(); }),
              expected, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.y() * o.y(); }),
              expected, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.z() * o.z(); }),
              expected, EXP_NEAR_TOLERANCE);
}

INSTANTIATE_TEST_SUITE_P(SufficientResolution, GLTMomentTest,
                         ::testing::Values(std::make_tuple(3, 2),
                                           std::make_tuple(3, 3),
                                           std::make_tuple(4, 2),
                                           std::make_tuple(5, 4),
                                           std::make_tuple(6, 6)));

// First and cross moments should vanish as soon as n_azim >= 2, well before
// n_azim is large enough for the diagonal second moments to be exact.
TEST(GLTMomentEdgeCaseTest, FirstAndCrossMomentsVanishEvenAtMinimalOrder) {
  GaussLegendreTrapezoid quad(2, 2);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.x(); }), 0.0,
              EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.y(); }), 0.0,
              EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.z(); }), 0.0,
              EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.x() * o.y(); }), 0.0,
              EXP_NEAR_TOLERANCE);
}

// At n_azim == 2 the trapezoid rule aliases the k=2 harmonic (n_azim divides
// k), so the x^2 and y^2 moments are individually wrong even though their
// sum (which only depends on the polar angle) is still correct. This
// documents that known limitation rather than treating it as a regression.
TEST(GLTMomentEdgeCaseTest,
     DiagonalMomentsAreAliasedWithOnlyTwoAzimuthalPoints) {
  GaussLegendreTrapezoid quad(2, 4);
  const double x_squared_moment =
      WeightedSum(quad, [](Ordinate o) { return o.x() * o.x(); });
  const double y_squared_moment =
      WeightedSum(quad, [](Ordinate o) { return o.y() * o.y(); });
  const double expected = kFourPi / 3.0;

  EXPECT_GT(std::abs(x_squared_moment - expected), 1e-3);
  EXPECT_GT(std::abs(y_squared_moment - expected), 1e-3);
  // The sum only depends on the polar angle (x^2 + y^2 = 1 - z^2), so it's
  // unaffected by azimuthal aliasing.
  EXPECT_NEAR(x_squared_moment + y_squared_moment, 2.0 * expected,
              EXP_NEAR_TOLERANCE);
}

// ---------------------------------------------------------------------------
// Geometric range checks
// ---------------------------------------------------------------------------

TEST(GLTRangeTest, zIsStrictlyInsideUnitInterval) {
  // Gauss-Legendre nodes never include the endpoints -1/1, so no ordinate's
  // z-direction cosine should be exactly +-1.
  GaussLegendreTrapezoid quad(4, 5);
  for (unsigned int i = 0; i < quad.n_points(); ++i) {
    Ordinate ord = quad.GetAbscissa(i);
    EXPECT_GT(ord.z(), -1.0);
    EXPECT_LT(ord.z(), 1.0);
  }
}

// ---------------------------------------------------------------------------
// Minimal/edge-case configuration
// ---------------------------------------------------------------------------

TEST(GLTEdgeCaseTest, SinglePointPerDimensionProducesOneOrdinate) {
  GaussLegendreTrapezoid quad(1, 1);
  EXPECT_EQ(quad.n_points(), 1u);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate) { return 1.0; }), kFourPi,
              EXP_NEAR_TOLERANCE);
}

// ---------------------------------------------------------------------------
// Exactness tests against closed-form sphere integrals
//
// These go beyond the moment checks above by verifying the quadrature
// against an independently-derived analytic formula for a whole family of
// monomials in x, y, z. Gauss-Legendre integrates a degree-d polynomial in
// mu exactly once 2*n_polar - 1 >= d, and the trapezoid rule integrates
// trigonometric polynomials of degree <= n_azim - 1 exactly, so with
// n_azim = 8 and n_polar = 6 every monomial up to total degree 6 below
// should be reproduced to quadrature tolerance.
// ---------------------------------------------------------------------------

class GLTMonomialExactnessTest : public ::testing::TestWithParam<Monomial> {};

TEST_P(GLTMonomialExactnessTest, MatchesClosedFormSphereIntegral) {
  const Monomial m = GetParam();
  const size_t n_azim = 8;
  const size_t n_polar = 6;
  GaussLegendreTrapezoid quad(n_azim, n_polar);

  const double expected = AnalyticSphereMonomialIntegral(m.a, m.b, m.c);
  const double actual = WeightedSum(quad, [&](Ordinate o) {
    return std::pow(o.x(), m.a) * std::pow(o.y(), m.b) * std::pow(o.z(), m.c);
  });
  EXPECT_NEAR(actual, expected, EXP_NEAR_TOLERANCE)
      << "Monomial x^" << m.a << " y^" << m.b << " z^" << m.c
      << " did not match its closed-form sphere integral.";
}

INSTANTIATE_TEST_SUITE_P(
    UpToDegreeSix, GLTMonomialExactnessTest,
    ::testing::Values(Monomial{0, 0, 0}, Monomial{2, 0, 0}, Monomial{0, 2, 0},
                      Monomial{0, 0, 2}, Monomial{2, 2, 0}, Monomial{2, 0, 2},
                      Monomial{0, 2, 2}, Monomial{4, 0, 0}, Monomial{0, 4, 0},
                      Monomial{0, 0, 4}, Monomial{4, 2, 0}, Monomial{2, 4, 0},
                      Monomial{0, 0, 6}, Monomial{2, 2, 2}));

// A monomial with azimuthal degree (a + b) at or beyond n_azim aliases,
// exactly like the diagonal-second-moment edge case above -- this checks
// that the exactness result genuinely depends on resolution rather than
// happening to hold for any inputs.
TEST(GLTMonomialExactnessTest, DegreeExceedingNAzimIsNotExact) {
  const size_t n_azim = 4;
  const size_t n_polar = 6;
  GaussLegendreTrapezoid quad(n_azim, n_polar);

  // x^4 has azimuthal content up to cos(4*phi), which aliases when
  // n_azim == 4.
  const double expected = AnalyticSphereMonomialIntegral(4, 0, 0);
  const double actual =
      WeightedSum(quad, [](Ordinate o) { return std::pow(o.x(), 4); });
  EXPECT_GT(std::abs(actual - expected), 1e-3);
}

// ---------------------------------------------------------------------------
// Convergence tests against non-polynomial integrands
//
// Monomial exactness only demonstrates correctness for functions the
// quadrature is designed to integrate exactly. These tests check the
// quadrature against smooth but non-polynomial integrands with known
// closed-form integrals, and confirm that increasing resolution actually
// drives the error down (rather than just checking a single fixed order).
// ---------------------------------------------------------------------------

TEST(GLTConvergenceTest, PolarErrorShrinksAsNPolarIncreasesForExpMu) {
  // int_{S^2} exp(mu) dOmega = 2*pi * (e - 1/e), independent of any
  // trigonometric-exactness concerns since it doesn't depend on phi at all.
  // Gauss-Legendre has no finite polynomial degree that reproduces exp(mu)
  // exactly, so this checks genuine convergence rather than exactness.
  const double expected = kTwoPi * (std::exp(1.0) - std::exp(-1.0));
  const std::vector<size_t> polar_orders = {2, 4, 6, 8};
  const size_t n_azim =
      2;  // Irrelevant here: the integrand is phi-independent.

  double previous_error = std::numeric_limits<double>::infinity();
  for (size_t n_polar : polar_orders) {
    GaussLegendreTrapezoid quad(n_azim, n_polar);
    const double actual =
        WeightedSum(quad, [](Ordinate o) { return std::exp(o.z()); });
    const double error = std::abs(actual - expected);
    EXPECT_LE(error, previous_error)
        << "Error did not shrink going to n_polar = " << n_polar;
    previous_error = error;
  }
  EXPECT_LT(previous_error, 1e-10)
      << "Expected near machine-precision agreement at the highest order.";
}

TEST(GLTConvergenceTest, AzimuthalErrorShrinksAsNAzimIncreasesForExpCosPhi) {
  // int_{S^2} exp(cos(phi)) dOmega = 4*pi*I_0(1), where I_0 is the modified
  // Bessel function of the first kind. exp(cos(phi)) is smooth and periodic
  // in phi, so the trapezoid rule should converge toward it very quickly as
  // n_azim grows, while n_polar plays no role since the integrand doesn't
  // depend on mu.
  constexpr double kBesselI0OfOne = 1.2660658777520084;
  const double expected = kFourPi * kBesselI0OfOne;
  const std::vector<size_t> azim_orders = {2, 4, 6, 8, 12};
  const size_t n_polar =
      3;  // Irrelevant here: the integrand is mu-independent.

  double previous_error = std::numeric_limits<double>::infinity();
  for (size_t n_azim : azim_orders) {
    GaussLegendreTrapezoid quad(n_azim, n_polar);
    const double actual = WeightedSum(
        quad, [](Ordinate o) { return std::exp(std::cos(Phi(o))); });
    const double error = std::abs(actual - expected);
    EXPECT_LE(error, previous_error)
        << "Error did not shrink going to n_azim = " << n_azim;
    previous_error = error;
  }
  EXPECT_LT(previous_error, 1e-10)
      << "Expected near machine-precision agreement at the highest order.";
}

}  // namespace hummingbird