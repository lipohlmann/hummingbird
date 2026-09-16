#include <gtest/gtest.h>

#include <cmath>

#include "quadrature/angular/ordinate.h"
#include "utils/constants.h"

namespace hummingbird {

// ---------------------------------------------------------------------------
// Ordinate::Reflect
// ---------------------------------------------------------------------------

TEST(OrdinateReflectTest, NormalIncidenceReversesDirection) {
  // Direction exactly along the surface normal bounces straight back.
  const arma::vec3 normal = {0.0, 1.0, 0.0};
  Ordinate incident(M_PI / 2.0, M_PI / 2.0);  // (0, 1, 0)
  ASSERT_NEAR(incident.x(), 0.0, EXP_NEAR_TOLERANCE);
  ASSERT_NEAR(incident.y(), 1.0, EXP_NEAR_TOLERANCE);
  ASSERT_NEAR(incident.z(), 0.0, EXP_NEAR_TOLERANCE);

  Ordinate reflected = Ordinate::Reflect(incident, normal);

  EXPECT_NEAR(reflected.x(), 0.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(reflected.y(), -1.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(reflected.z(), 0.0, EXP_NEAR_TOLERANCE);
}

TEST(OrdinateReflectTest, GrazingIncidenceIsUnaffected) {
  // Direction perpendicular to the normal (d.n == 0) passes through
  // unchanged.
  const arma::vec3 normal = {0.0, 0.0, 1.0};
  Ordinate incident(0.0, M_PI / 2.0);  // (1, 0, 0)
  ASSERT_NEAR(incident.x(), 1.0, EXP_NEAR_TOLERANCE);
  ASSERT_NEAR(incident.y(), 0.0, EXP_NEAR_TOLERANCE);
  ASSERT_NEAR(incident.z(), 0.0, EXP_NEAR_TOLERANCE);

  Ordinate reflected = Ordinate::Reflect(incident, normal);

  EXPECT_NEAR(reflected.x(), 1.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(reflected.y(), 0.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(reflected.z(), 0.0, EXP_NEAR_TOLERANCE);
}

TEST(OrdinateReflectTest, FlipsXDirectionCosineAcross1DEastBoundary) {
  // The case this exists for: a 1D boundary with outward normal (1,0,0)
  // (e.g. "bc:east") reflects an outgoing mu=1 direction straight back as
  // mu=-1.
  const arma::vec3 normal = {1.0, 0.0, 0.0};
  Ordinate incident(0.0, M_PI / 2.0);  // (1, 0, 0)

  Ordinate reflected = Ordinate::Reflect(incident, normal);

  EXPECT_NEAR(reflected.x(), -1.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(reflected.y(), 0.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(reflected.z(), 0.0, EXP_NEAR_TOLERANCE);
}

TEST(OrdinateReflectTest, ObliqueIncidenceMatchesHandDerivedDirection) {
  // 45 degrees off the normal in the x-z plane: reflecting across the
  // z-normal should flip only the z-component.
  const arma::vec3 normal = {0.0, 0.0, 1.0};
  const double s = std::sqrt(2.0) / 2.0;
  Ordinate incident(0.0, M_PI / 4.0);  // (s, 0, s)
  ASSERT_NEAR(incident.x(), s, EXP_NEAR_TOLERANCE);
  ASSERT_NEAR(incident.z(), s, EXP_NEAR_TOLERANCE);

  Ordinate reflected = Ordinate::Reflect(incident, normal);

  EXPECT_NEAR(reflected.x(), s, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(reflected.y(), 0.0, EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(reflected.z(), -s, EXP_NEAR_TOLERANCE);
}

TEST(OrdinateReflectTest, ReflectingTwiceAcrossSameNormalRestoresOriginal) {
  const arma::vec3 normal = arma::normalise(arma::vec3{1.0, 1.0, 1.0});
  Ordinate original(0.3, 1.2);

  Ordinate once = Ordinate::Reflect(original, normal);
  Ordinate twice = Ordinate::Reflect(once, normal);

  EXPECT_NEAR(twice.x(), original.x(), EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(twice.y(), original.y(), EXP_NEAR_TOLERANCE);
  EXPECT_NEAR(twice.z(), original.z(), EXP_NEAR_TOLERANCE);
}

TEST(OrdinateReflectTest, ResultIsAUnitVector) {
  const arma::vec3 normal = arma::normalise(arma::vec3{0.2, -0.5, 0.8});
  Ordinate incident(-0.7, 2.1);

  Ordinate reflected = Ordinate::Reflect(incident, normal);

  double norm_sq = reflected.x() * reflected.x() +
                   reflected.y() * reflected.y() +
                   reflected.z() * reflected.z();
  EXPECT_NEAR(norm_sq, 1.0, EXP_NEAR_TOLERANCE);
}

}  // namespace hummingbird
