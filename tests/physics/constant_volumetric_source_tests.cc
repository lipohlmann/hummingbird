#include <gtest/gtest.h>

#include <cmath>

#include "physics/constant_volumetric_source.h"

namespace hummingbird {
TEST(ConstantVolumetricSourceTest, NoInputNeeded) {
  ConstantVolumetricSource cvs(4.0 * M_PI);
  EXPECT_DOUBLE_EQ(cvs.EvaluateAtNode(), 1.0);
}

TEST(ConstantVolumetricSourceTest, InputIgnored) {
  Node dummy_node(1, 5, 2, 4);
  Ordinate dummy_ordinate(M_PI / 4.0, M_PI / 4.0);
  ConstantVolumetricSource cvs(4.0 * M_PI);
  EXPECT_DOUBLE_EQ(cvs.EvaluateAtNode(dummy_node), 1.0);
  EXPECT_DOUBLE_EQ(cvs.EvaluateAtNode(dummy_node, dummy_ordinate), 1.0);
}
}  // namespace hummingbird
