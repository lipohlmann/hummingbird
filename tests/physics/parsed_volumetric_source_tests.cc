#include <gtest/gtest.h>

#include <string>

#include "mesh/node.h"
#include "physics/parsed_volumetric_source.h"
#include "quadrature/angular/ordinate.h"

namespace hummingbird {
class ParsedVolumetricSourceTest : public testing::Test {
 protected:
  Node node{.id = 1, .x = 1, .y = 2, .z = 5};
};

TEST_F(ParsedVolumetricSourceTest, EvaluateAtNode) {
  Ordinate ordinate(0, 0);
  std::string expression = "x+y+z";
  ParsedVolumetricSource pvs(expression);
  EXPECT_DOUBLE_EQ(pvs.EvaluateAtNode(node, ordinate), 8);
}
}  // namespace hummingbird
