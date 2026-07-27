#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <tuple>
#include <vector>

#include "quadrature/angular/gauss_legendre_chebyshev.h"
#include "quadrature/angular/ordinate.h"
#include "quadrature/gauss_chebyshev.h"
#include "quadrature/gauss_legendre.h"
#include "utils/constants.h"

namespace hummingbird::quadrature::angular {

namespace {

constexpr double kFourPi = 4.0 * M_PI;

// Sums a per-ordinate quantity over every point in the quadrature, weighted
// by that point's quadrature weight. This plays the same role as
// QuadratureBase::Integrate(), but works directly off direction cosines
// since Ordinate (unlike double) isn't a scalar abscissa type the shared
// helper machinery expects.
double WeightedSum(const GaussLegendreChebyshev& quad,
                   const std::function<double(Ordinate)>& f) {
  double sum = 0.0;
  for (unsigned int i = 0; i < quad.n_points(); ++i) {
    Ordinate ord = quad.GetAbscissa(i);
    sum += quad.GetWeight(i) * f(ord);
  }
  return sum;
}

}  // namespace

// ---------------------------------------------------------------------------
// Structural tests
// ---------------------------------------------------------------------------

class GLCStructureTest
    : public ::testing::TestWithParam<std::tuple<size_t, size_t>> {};

TEST_P(GLCStructureTest, ProducesExpectedNumberOfOrdinates) {
  auto [n_half_azim, n_polar] = GetParam();
  GaussLegendreChebyshev quad(n_half_azim, n_polar);
  EXPECT_EQ(quad.n_points(), 2 * n_half_azim * n_polar);
}

TEST_P(GLCStructureTest, EveryOrdinateIsAUnitVector) {
  auto [n_half_azim, n_polar] = GetParam();
  GaussLegendreChebyshev quad(n_half_azim, n_polar);
  for (unsigned int i = 0; i < quad.n_points(); ++i) {
    Ordinate ord = quad.GetAbscissa(i);
    const double norm_sq = ord.XCosine() * ord.XCosine() +
                           ord.YCosine() * ord.YCosine() +
                           ord.ZCosine() * ord.ZCosine();
    EXPECT_NEAR(norm_sq, 1.0, utils::EXP_NEAR_TOLERANCE)
        << "Direction cosines at index " << i << " are not normalized.";
  }
}

TEST_P(GLCStructureTest, AllWeightsArePositive) {
  auto [n_half_azim, n_polar] = GetParam();
  GaussLegendreChebyshev quad(n_half_azim, n_polar);
  for (unsigned int i = 0; i < quad.n_points(); ++i) {
    EXPECT_GT(quad.GetWeight(i), 0.0) << "Non-positive weight at index " << i;
  }
}

TEST_P(GLCStructureTest, WeightsSumToFourPi) {
  // The quadrature approximates the surface integral over the unit sphere,
  // whose total solid angle is 4*pi. This holds for any n_half_azim/n_polar
  // combination, since it only relies on the Gauss-Legendre weights summing
  // to 2 (length of [-1,1]) and the mirrored Gauss-Chebyshev weights summing
  // to 2*pi (full azimuthal circle).
  auto [n_half_azim, n_polar] = GetParam();
  GaussLegendreChebyshev quad(n_half_azim, n_polar);
  const double sum = WeightedSum(quad, [](Ordinate) { return 1.0; });
  EXPECT_NEAR(sum, kFourPi, utils::EXP_NEAR_TOLERANCE * 10.0);
}

TEST_P(GLCStructureTest, IntegrateAgreesWithManualWeightedSumForAConstant) {
  // Exercise the inherited QuadratureBase::Integrate() / QuadraturePair path
  // directly, since that's the interface production code will actually use.
  auto [n_half_azim, n_polar] = GetParam();
  GaussLegendreChebyshev quad(n_half_azim, n_polar);

  std::vector<QuadraturePair> pairs;
  pairs.reserve(quad.n_points());
  for (unsigned int i = 0; i < quad.n_points(); ++i) {
    pairs.push_back({static_cast<int>(i), 1.0});
  }

  EXPECT_NEAR(quad.Integrate(pairs), kFourPi, utils::EXP_NEAR_TOLERANCE);
}

TEST_P(GLCStructureTest, WeightsMatchProductOfComponentQuadratures) {
  // Each ordinate's weight should be w_polar(i) * w_azim(j), independent of
  // whether it's the "positive" or mirrored "negative" azimuthal copy.
  auto [n_half_azim, n_polar] = GetParam();
  GaussLegendreChebyshev quad(n_half_azim, n_polar);
  GaussChebyshev gc_azim(n_half_azim);
  GaussLegendre gl_polar(n_polar);

  unsigned int counter = 0;
  for (unsigned int i = 0; i < n_polar; ++i) {
    for (unsigned int j = 0; j < n_half_azim; ++j) {
      const double expected = gl_polar.GetWeight(i) * gc_azim.GetWeight(j);
      EXPECT_NEAR(quad.GetWeight(counter), expected, utils::EXP_NEAR_TOLERANCE)
          << "Mismatch at positive ordinate i=" << i << ", j=" << j;
      ++counter;
      EXPECT_NEAR(quad.GetWeight(counter), expected, utils::EXP_NEAR_TOLERANCE)
          << "Mismatch at mirrored ordinate i=" << i << ", j=" << j;
      ++counter;
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    VariousOrders, GLCStructureTest,
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
// aliasing (n_half_azim >= 2 for the azimuthal moments, since the mirrored
// Gauss-Chebyshev construction reduces to a 2*n_half_azim-point equally
// spaced rule on the circle, and n_polar >= 2 for the polar moment, since
// mu^2 is a degree-2 polynomial).
// ---------------------------------------------------------------------------

class GLCMomentTest
    : public ::testing::TestWithParam<std::tuple<size_t, size_t>> {};

TEST_P(GLCMomentTest, FirstMomentsVanish) {
  auto [n_half_azim, n_polar] = GetParam();
  GaussLegendreChebyshev quad(n_half_azim, n_polar);

  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.XCosine(); }), 0.0,
              utils::EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.YCosine(); }), 0.0,
              utils::EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate o) { return o.ZCosine(); }), 0.0,
              utils::EXP_NEAR_TOLERANCE);
}

TEST_P(GLCMomentTest, CrossMomentsVanish) {
  auto [n_half_azim, n_polar] = GetParam();
  GaussLegendreChebyshev quad(n_half_azim, n_polar);

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

TEST_P(GLCMomentTest, DiagonalSecondMomentsEqualFourPiOverThree) {
  auto [n_half_azim, n_polar] = GetParam();
  GaussLegendreChebyshev quad(n_half_azim, n_polar);
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

INSTANTIATE_TEST_SUITE_P(SufficientResolution, GLCMomentTest,
                         ::testing::Values(std::make_tuple(2, 2),
                                           std::make_tuple(3, 3),
                                           std::make_tuple(4, 2),
                                           std::make_tuple(2, 4),
                                           std::make_tuple(6, 6)));

// First and cross moments should vanish even when resolution is too coarse
// for the diagonal second moments to be exact, since that cancellation comes
// from quadrature symmetry, not point count.
TEST(GLCMomentEdgeCaseTest, FirstAndCrossMomentsVanishEvenAtMinimalOrder) {
  GaussLegendreChebyshev quad(1, 2);
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

// At n_half_azim == 1 there are only two azimuthal directions (phi = pi/2
// and phi = 3*pi/2), so cos(phi) == 0 identically and the x^2 moment is
// aliased to zero rather than the true 4*pi/3. This documents that known
// limitation rather than treating it as a regression.
TEST(GLCMomentEdgeCaseTest, XSquaredMomentIsAliasedWithOnlyOneAzimuthalPoint) {
  GaussLegendreChebyshev quad(1, 4);
  const double x_squared_moment =
      WeightedSum(quad, [](Ordinate o) { return o.XCosine() * o.XCosine(); });
  EXPECT_NEAR(x_squared_moment, 0.0, utils::EXP_NEAR_TOLERANCE);
  EXPECT_GT(std::abs(x_squared_moment - kFourPi / 3.0), 1e-3);
}

// ---------------------------------------------------------------------------
// Geometric range checks
// ---------------------------------------------------------------------------

TEST(GLCRangeTest, ZCosineIsStrictlyInsideUnitInterval) {
  // Gauss-Legendre nodes never include the endpoints -1/1, so no ordinate's
  // z-direction cosine should be exactly +-1.
  GaussLegendreChebyshev quad(4, 5);
  for (unsigned int i = 0; i < quad.n_points(); ++i) {
    Ordinate ord = quad.GetAbscissa(i);
    EXPECT_GT(ord.ZCosine(), -1.0);
    EXPECT_LT(ord.ZCosine(), 1.0);
  }
}

TEST(GLCRangeTest, MirroredPairsShareThePolarAngleButOpposeAzimuthally) {
  // Ordinates are pushed in (positive, negative) pairs for each (i, j): the
  // pair shares a z-cosine (polar angle) and an x-cosine (cos(phi) ==
  // cos(2*pi - phi)), but differ in sign on the y-cosine, since only the
  // azimuthal angle is mirrored (phi -> 2*pi - phi).
  const size_t n_half_azim = 3;
  const size_t n_polar = 2;
  GaussLegendreChebyshev quad(n_half_azim, n_polar);

  for (unsigned int idx = 0; idx < quad.n_points(); idx += 2) {
    Ordinate positive = quad.GetAbscissa(idx);
    Ordinate negative = quad.GetAbscissa(idx + 1);

    EXPECT_NEAR(positive.ZCosine(), negative.ZCosine(),
                utils::EXP_NEAR_TOLERANCE)
        << "Mirrored pair at index " << idx << " should share a polar angle.";
    EXPECT_NEAR(positive.XCosine(), negative.XCosine(),
                utils::EXP_NEAR_TOLERANCE)
        << "Mirrored pair at index " << idx
        << " should share an x-direction cosine.";
    EXPECT_NEAR(positive.YCosine(), -negative.YCosine(),
                utils::EXP_NEAR_TOLERANCE)
        << "Mirrored pair at index " << idx
        << " should have opposite y-direction cosines.";
    EXPECT_NEAR(quad.GetWeight(idx), quad.GetWeight(idx + 1),
                utils::EXP_NEAR_TOLERANCE)
        << "Mirrored pair at index " << idx << " should share a weight.";
  }
}

// ---------------------------------------------------------------------------
// Minimal/edge-case configuration
// ---------------------------------------------------------------------------

TEST(GLCEdgeCaseTest, SinglePointPerHalfProducesTwoOrdinates) {
  GaussLegendreChebyshev quad(1, 1);
  EXPECT_EQ(quad.n_points(), 2u);
  EXPECT_NEAR(WeightedSum(quad, [](Ordinate) { return 1.0; }), kFourPi,
              utils::EXP_NEAR_TOLERANCE);
}

}  // namespace hummingbird::quadrature::angular