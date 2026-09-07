#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include "mesh/node.h"
#include "physics/parsed_volumetric_source.h"
#include "quadrature/angular/ordinate.h"

namespace hummingbird {

// NOTE: ParsedVolumetricSource::EvaluateAtNode uses mu = ordinate.x(), which
// is a direction cosine, not a free parameter. Per Ordinate's documented
// convention (azimuth in [-pi,pi], polar in [0,pi]):
//   x() = sin(polar) * cos(azimuth)
//   y() = sin(polar) * sin(azimuth)
//   z() = cos(polar)
// Test angles below are chosen to land on exact/near-exact x() values
// (0, +-0.5, +-1) rather than round-tripping through asin/sin.

class ParsedVolumetricSourceTest : public testing::Test {
 protected:
  Node node{.id = 1, .x = 1, .y = 2, .z = 5};
};

// ---------------------------------------------------------------------------
// Baseline / existing behavior
// ---------------------------------------------------------------------------

TEST_F(ParsedVolumetricSourceTest, EvaluateAtNode) {
  Ordinate ordinate(0, 0);  // mu unused by this expression.
  std::string expression = "x+y+z";
  ParsedVolumetricSource pvs(expression);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), 8);
}

// ---------------------------------------------------------------------------
// Individual variable wiring
// ---------------------------------------------------------------------------

TEST_F(ParsedVolumetricSourceTest, ReadsXOnly) {
  ParsedVolumetricSource pvs("x");
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), node.x);
}

TEST_F(ParsedVolumetricSourceTest, ReadsYOnly) {
  ParsedVolumetricSource pvs("y");
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), node.y);
}

TEST_F(ParsedVolumetricSourceTest, ReadsZOnly) {
  ParsedVolumetricSource pvs("z");
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), node.z);
}

TEST_F(ParsedVolumetricSourceTest, ReadsMuAtEquatorZeroAzimuth) {
  // polar = pi/2, azimuth = 0 => x() = sin(pi/2)*cos(0) = 1
  ParsedVolumetricSource pvs("mu");
  Ordinate ordinate(0, M_PI / 2.0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), 1.0);
}

TEST_F(ParsedVolumetricSourceTest, ReadsMuAtPoleIsZeroRegardlessOfAzimuth) {
  // polar = 0 => x() = sin(0)*cos(azimuth) = 0, for any azimuth.
  ParsedVolumetricSource pvs("mu");
  Ordinate ordinate(1.2345, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), 0.0);
}

TEST_F(ParsedVolumetricSourceTest, ReadsPositiveMu) {
  // polar = pi/6, azimuth = 0 => x() = sin(pi/6)*cos(0) = 0.5
  ParsedVolumetricSource pvs("mu");
  Ordinate ordinate(0, M_PI / 6.0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), 0.5);
}

TEST_F(ParsedVolumetricSourceTest, ReadsNegativeMu) {
  // polar = pi/6, azimuth = pi => x() = sin(pi/6)*cos(pi) = -0.5
  ParsedVolumetricSource pvs("mu");
  Ordinate ordinate(M_PI, M_PI / 6.0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), -0.5);
}

// ---------------------------------------------------------------------------
// Coordinates with negative / zero values
// ---------------------------------------------------------------------------

TEST_F(ParsedVolumetricSourceTest, HandlesNegativeCoordinates) {
  Node n{.id = 2, .x = -3, .y = -4, .z = -5};
  ParsedVolumetricSource pvs("x+y+z");
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(n, ordinate), -12);
}

TEST_F(ParsedVolumetricSourceTest, HandlesZeroCoordinatesAndZeroMu) {
  Node n{.id = 3, .x = 0, .y = 0, .z = 0};
  ParsedVolumetricSource pvs("x+y+z+mu");
  Ordinate ordinate(0, 0);  // mu = 0
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(n, ordinate), 0);
}

// ---------------------------------------------------------------------------
// Operators and precedence (mu unused; Ordinate arbitrary but valid)
// ---------------------------------------------------------------------------

TEST_F(ParsedVolumetricSourceTest, Multiplication) {
  ParsedVolumetricSource pvs("x*y*z");
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), 10);  // 1*2*5
}

TEST_F(ParsedVolumetricSourceTest, Division) {
  ParsedVolumetricSource pvs("z/y");
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), 2.5);
}

TEST_F(ParsedVolumetricSourceTest, Exponentiation) {
  ParsedVolumetricSource pvs("y^3");
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), 8);
}

TEST_F(ParsedVolumetricSourceTest, OperatorPrecedenceRespected) {
  ParsedVolumetricSource pvs("x+y*z");  // 1 + 10 = 11, not (1+2)*5 = 15
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), 11);
}

TEST_F(ParsedVolumetricSourceTest, ParenthesesOverridePrecedence) {
  ParsedVolumetricSource pvs("(x+y)*z");
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), 15);
}

TEST_F(ParsedVolumetricSourceTest, UnaryMinus) {
  ParsedVolumetricSource pvs("-x");
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), -1);
}

// ---------------------------------------------------------------------------
// Built-in math functions / constants
// ---------------------------------------------------------------------------

TEST_F(ParsedVolumetricSourceTest, SqrtFunction) {
  ParsedVolumetricSource pvs("sqrt(z)");
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), std::sqrt(5.0));
}

TEST_F(ParsedVolumetricSourceTest, TrigFunctionOfMu) {
  // polar = pi/6, azimuth = 0 => mu = 0.5
  ParsedVolumetricSource pvs("cos(mu)");
  Ordinate ordinate(0, M_PI / 6.0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), std::cos(0.5));
}

TEST_F(ParsedVolumetricSourceTest, ExpFunction) {
  ParsedVolumetricSource pvs("exp(x)");
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), std::exp(1.0));
}

TEST_F(ParsedVolumetricSourceTest, AbsFunctionOnNegativeInput) {
  Node n{.id = 4, .x = -7, .y = 0, .z = 0};
  ParsedVolumetricSource pvs("abs(x)");
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(n, ordinate), 7);
}

TEST_F(ParsedVolumetricSourceTest, PiConstantTimesMu) {
  // polar = pi/2, azimuth = 0 => mu = 1
  ParsedVolumetricSource pvs("pi*mu");
  Ordinate ordinate(0, M_PI / 2.0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), M_PI);
}

// ---------------------------------------------------------------------------
// Combined spatial + angular MMS-style expression
// ---------------------------------------------------------------------------

TEST_F(ParsedVolumetricSourceTest, CombinedSpatialAngularExpression) {
  // polar = pi/6, azimuth = 0 => mu = 0.5
  ParsedVolumetricSource pvs("sin(x)*cos(y)+mu*z");
  Ordinate ordinate(0, M_PI / 6.0);
  double expected = std::sin(1.0) * std::cos(2.0) + 0.5 * 5.0;
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), expected);
}

// ---------------------------------------------------------------------------
// Formatting robustness
// ---------------------------------------------------------------------------

TEST_F(ParsedVolumetricSourceTest, ToleratesWhitespace) {
  ParsedVolumetricSource pvs("  x  +  y  +  z  ");
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), 8);
}

TEST_F(ParsedVolumetricSourceTest, ConstantExpressionIgnoresNodeAndOrdinate) {
  ParsedVolumetricSource pvs("42");
  Ordinate ordinate(0, M_PI / 2.0);  // mu = 1, irrelevant to a literal
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), 42);
}

TEST_F(ParsedVolumetricSourceTest, ScientificNotationLiteral) {
  ParsedVolumetricSource pvs("1e-3*x");
  Ordinate ordinate(0, 0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), 1e-3);
}

// ---------------------------------------------------------------------------
// Reuse across multiple nodes/ordinates: confirms no state leaks between
// calls (a fresh te_parser + local bindings are built each call).
// ---------------------------------------------------------------------------

TEST_F(ParsedVolumetricSourceTest, RepeatedCallsAreIndependent) {
  ParsedVolumetricSource pvs("x+mu");
  Node n1{.id = 5, .x = 1, .y = 0, .z = 0};
  Node n2{.id = 6, .x = 100, .y = 0, .z = 0};
  Ordinate o1(0, 0);           // mu = 0
  Ordinate o2(0, M_PI / 2.0);  // mu = 1

  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(n1, o1), 1.0);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(n2, o2), 101.0);
  // Re-run the first pair to ensure no residual state from the second call.
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(n1, o1), 1.0);
}

// ---------------------------------------------------------------------------
// Error handling
// ---------------------------------------------------------------------------

TEST_F(ParsedVolumetricSourceTest, ThrowsOnMalformedExpression) {
  ParsedVolumetricSource pvs("x+*y");
  Ordinate ordinate(0, 0);
  EXPECT_THROW(pvs.EvaluateAtNode(node, ordinate), std::runtime_error);
}

TEST_F(ParsedVolumetricSourceTest, ThrowsOnUnbalancedParentheses) {
  ParsedVolumetricSource pvs("(x+y");
  Ordinate ordinate(0, 0);
  EXPECT_THROW(pvs.EvaluateAtNode(node, ordinate), std::runtime_error);
}

TEST_F(ParsedVolumetricSourceTest, ThrowsOnUndefinedVariable) {
  ParsedVolumetricSource pvs("x+t");  // "t" is not registered
  Ordinate ordinate(0, 0);
  EXPECT_THROW(pvs.EvaluateAtNode(node, ordinate), std::runtime_error);
}

TEST_F(ParsedVolumetricSourceTest, ThrowsOnEmptyExpression) {
  ParsedVolumetricSource pvs("");
  Ordinate ordinate(0, 0);
  EXPECT_THROW(pvs.EvaluateAtNode(node, ordinate), std::runtime_error);
}

TEST_F(ParsedVolumetricSourceTest, ThrowsOnUnknownFunctionName) {
  ParsedVolumetricSource pvs("frobnicate(x)");
  Ordinate ordinate(0, 0);
  EXPECT_THROW(pvs.EvaluateAtNode(node, ordinate), std::runtime_error);
}

TEST_F(ParsedVolumetricSourceTest, DivisionByZeroProducesInfNotThrow) {
  Node n{.id = 7, .x = 1, .y = 0, .z = 0};
  ParsedVolumetricSource pvs("x/y");
  Ordinate ordinate(0, 0);
  EXPECT_THROW(pvs.EvaluateAtNode(n, ordinate), std::runtime_error);
}

TEST_F(ParsedVolumetricSourceTest, ZeroOverZeroProducesNan) {
  Node n{.id = 8, .x = 0, .y = 0, .z = 0};
  ParsedVolumetricSource pvs("x/y");
  Ordinate ordinate(0, 0);
  EXPECT_THROW(pvs.EvaluateAtNode(n, ordinate), std::runtime_error);
}

}  // namespace hummingbird