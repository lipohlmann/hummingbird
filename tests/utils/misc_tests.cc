#include <gtest/gtest.h>

#include "utils/misc.h"

namespace hummingbird::utils {
TEST(MathUtils, DoubleEqual) {
  EXPECT_TRUE(DoubleEqual(1.0, 1.0));
  EXPECT_TRUE(DoubleEqual(1.0, 1.0 + TOLERANCE / 2.0));
  EXPECT_FALSE(DoubleEqual(1.0, -1.0));
  EXPECT_FALSE(DoubleEqual(1.0, 2.5));
  EXPECT_FALSE(DoubleEqual(1.0, 1.0 + TOLERANCE));
}
}  // namespace hummingbird::utils
