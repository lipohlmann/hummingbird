#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <tuple>
#include <vector>

#include "quadrature/angular/gauss_legendre_trapezoid.h"
#include "quadrature/angular/ordinate.h"
#include "quadrature/gauss_legendre.h"
#include "utils/constants.h"

namespace hummingbird::quadrature::angular {

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

// Groups ordinate indices by shared z-cosine (i.e. by polar level), without
// assuming anything about how the implementation orders its points. Indices
// within a group are all the azimuthal samples at that polar angle.
std::vector<std::vector<size_t>> GroupIndicesByZCosine(
    const GaussLegendreTrapezoid& quad) {
  std::vector<size_t> order(quad.n_points());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    return quad.GetAbscissa(a).ZCosine() < quad.GetAbscissa(b).ZCosine();
  });

  std::vector<std::vector<size_t>> groups;
  for (size_t idx : order) {
    const double z = quad.GetAbscissa(idx).ZCosine();
    if (groups.empty() ||
        std::abs(quad.GetAbscissa(groups.back().front()).ZCosine() - z) >
            1e-6) {
      groups.push_back({idx});
    } else {
      groups.back().push_back(idx);
    }
  }
  return groups;
}

// Azimuthal angle in [0, 2*pi).
double Phi(Ordinate ord) {
  double phi = std::atan2(ord.YCosine(), ord.XCosine());
  if (phi < 0.0) phi += kTwoPi;
  return phi;
}

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
    const double norm_sq = ord.XCosine() * ord.XCosine() +
                           ord.YCosine() * ord.YCosine() +
                           ord.ZCosine() * ord.ZCosine();
    EXPECT_NEAR(norm_sq, 1.0, utils::EXP_NEAR_TOLERANCE)
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
  EXPECT_NEAR(sum, kFourPi, utils::EXP_NEAR_TOLERANCE * 10.0);
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

  EXPECT_NEAR(quad.Integrate(pairs), kFourPi, utils::EXP_NEAR_TOLERANCE);
}

TEST_P(GLTStructureTest,
       WeightsMatchProductOfPolarWeightAndUniformAzimuthalWeight) {
  // Each ordinate's weight should be w_polar(i) * (2*pi / n_azim), since the
  // trapezoid rule on a full period assigns every azimuthal sample an equal
  // weight. Grouping by shared z-cosine avoids assuming anything about the
  // implementation's point ordering.
  auto [n_azim, n_polar] = GetParam();
  GaussLegendreTrapezoid quad(n_azim, n_polar);
  GaussLegendre gl_polar(n_polar);

  auto groups = GroupIndicesByZCosine(quad);
  ASSERT_EQ(groups.size(), n_polar);

  std::vector<double> actual_weights;
  for (const auto& group : groups) {
    ASSERT_EQ(group.size(), n_azim);
    const double first_weight = quad.GetWeight(group.front());
    for (size_t idx : group) {
      EXPECT_NEAR(quad.GetWeight(idx), first_weight, utils::EXP_NEAR_TOLERANCE)
          << "Azimuthal weights at a fixed polar level should be equal.";
    }
    actual_weights.push_back(first_weight);
  }
  std::sort(actual_weights.begin(), actual_weights.end());

  std::vector<double> expected_weights;
  for (unsigned int i = 0; i < n_polar; ++i) {
    expected_weights.push_back(gl_polar.GetWeight(i) * (kTwoPi / n_azim));
  }
  std::sort(expected_weights.begin(), expected_weights.end());

  ASSERT_EQ(actual_weights.size(), expected_weights.size());
  for (size_t i = 0; i < actual_weights.size(); ++i) {
    EXPECT_NEAR(actual_weights.at(i), expected_weights.at(i),
                utils::EXP_NEAR_TOLERANCE)
        << "Mismatch at sorted polar-weight index " << i;
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

  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.XCosine(); }), 0.0,
              utils::EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.YCosine(); }), 0.0,
              utils::EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.ZCosine(); }), 0.0,
              utils::EXP_NEAR_TOLERANCE);
}

TEST_P(GLTMomentTest, CrossMomentsVanish) {
  auto [n_azim, n_polar] = GetParam();
  GaussLegendreTrapezoid quad(n_azim, n_polar);

  EXPECT_NEAR(
      WeightedSum(quad, [](Ordinate o) { return o.XCosine() * o.YCosine(); }),
      0.0, utils::EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(
      WeightedSum(quad, [](Ordinate o) { return o.YCosine() * o.ZCosine(); }),
      0.0, utils::EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(
      WeightedSum(quad, [](Ordinate o) { return o.XCosine() * o.ZCosine(); }),
      0.0, utils::EXP_NEAR_TOLERANCE);
}

TEST_P(GLTMomentTest, DiagonalSecondMomentsEqualFourPiOverThree) {
  auto [n_azim, n_polar] = GetParam();
  GaussLegendreTrapezoid quad(n_azim, n_polar);
  const double expected = kFourPi / 3.0;

  EXPECT_NEAR(
      WeightedSum(quad, [](Ordinate o) { return o.XCosine() * o.XCosine(); }),
      expected, utils::EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(
      WeightedSum(quad, [](Ordinate o) { return o.YCosine() * o.YCosine(); }),
      expected, utils::EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(
      WeightedSum(quad, [](Ordinate o) { return o.ZCosine() * o.ZCosine(); }),
      expected, utils::EXP_NEAR_TOLERANCE);
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
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.XCosine(); }), 0.0,
              utils::EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.YCosine(); }), 0.0,
              utils::EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.ZCosine(); }), 0.0,
              utils::EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(
      WeightedSum(quad, [](Ordinate o) { return o.XCosine() * o.YCosine(); }),
      0.0, utils::EXP_NEAR_TOLERANCE);
}

// At n_azim == 2 the trapezoid rule aliases the k=2 harmonic (n_azim divides
// k), so the x^2 and y^2 moments are individually wrong even though their
// sum (which only depends on the polar angle) is still correct. This
// documents that known limitation rather than treating it as a regression.
TEST(GLTMomentEdgeCaseTest,
     DiagonalMomentsAreAliasedWithOnlyTwoAzimuthalPoints) {
  GaussLegendreTrapezoid quad(2, 4);
  const double x_squared_moment =
      WeightedSum(quad, [](Ordinate o) { return o.XCosine() * o.XCosine(); });
  const double y_squared_moment =
      WeightedSum(quad, [](Ordinate o) { return o.YCosine() * o.YCosine(); });
  const double expected = kFourPi / 3.0;

  EXPECT_GT(std::abs(x_squared_moment - expected), 1e-3);
  EXPECT_GT(std::abs(y_squared_moment - expected), 1e-3);
  // The sum only depends on the polar angle (x^2 + y^2 = 1 - z^2), so it's
  // unaffected by azimuthal aliasing.
  EXPECT_NEAR(x_squared_moment + y_squared_moment, 2.0 * expected,
              utils::EXP_NEAR_TOLERANCE);
}

// ---------------------------------------------------------------------------
// Geometric range checks
// ---------------------------------------------------------------------------

TEST(GLTRangeTest, ZCosineIsStrictlyInsideUnitInterval) {
  // Gauss-Legendre nodes never include the endpoints -1/1, so no ordinate's
  // z-direction cosine should be exactly +-1.
  GaussLegendreTrapezoid quad(4, 5);
  for (unsigned int i = 0; i < quad.n_points(); ++i) {
    Ordinate ord = quad.GetAbscissa(i);
    EXPECT_GT(ord.ZCosine(), -1.0);
    EXPECT_LT(ord.ZCosine(), 1.0);
  }
}

TEST(GLTRangeTest, PointsAtFixedPolarLevelAreEquallySpacedInAzimuth) {
  // For each polar level, the n_azim trapezoid samples should share a
  // z-cosine and be spaced 2*pi/n_azim apart in azimuthal angle. Points are
  // grouped by z-cosine rather than by assumed index layout, so this holds
  // regardless of how the implementation orders its points.
  const size_t n_azim = 4;
  const size_t n_polar = 3;
  GaussLegendreTrapezoid quad(n_azim, n_polar);

  auto groups = GroupIndicesByZCosine(quad);
  ASSERT_EQ(groups.size(), n_polar);

  const double expected_spacing = kTwoPi / n_azim;
  for (const auto& group : groups) {
    ASSERT_EQ(group.size(), n_azim);

    std::vector<double> phis;
    for (size_t idx : group) phis.push_back(Phi(quad.GetAbscissa(idx)));
    std::sort(phis.begin(), phis.end());

    for (size_t k = 0; k + 1 < phis.size(); ++k) {
      EXPECT_NEAR(phis.at(k + 1) - phis.at(k), expected_spacing,
                  utils::EXP_NEAR_TOLERANCE)
          << "Uneven azimuthal spacing within a polar level.";
    }
    // Wraparound gap between the last and first (plus a full turn).
    EXPECT_NEAR(phis.front() + kTwoPi - phis.back(), expected_spacing,
                utils::EXP_NEAR_TOLERANCE)
        << "Uneven azimuthal spacing across the 0/2*pi wraparound.";
  }
}

// ---------------------------------------------------------------------------
// Minimal/edge-case configuration
// ---------------------------------------------------------------------------

TEST(GLTEdgeCaseTest, SinglePointPerDimensionProducesOneOrdinate) {
  GaussLegendreTrapezoid quad(1, 1);
  EXPECT_EQ(quad.n_points(), 1u);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate) { return 1.0; }), kFourPi,
              utils::EXP_NEAR_TOLERANCE);
}

}  // namespace hummingbird::quadrature::angular