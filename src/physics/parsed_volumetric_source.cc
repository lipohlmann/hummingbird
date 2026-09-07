#include "physics/parsed_volumetric_source.h"

#include <stdexcept>

#include "utils/tinyexpr.h"

namespace hummingbird {
ParsedVolumetricSource::ParsedVolumetricSource(
    const std::string_view expression)
    : expression_(expression) {}

double ParsedVolumetricSource::EvaluateAtNode(const Node& node,
                                              const Ordinate& ordinate) {
  te_parser tep;
  te_type x{node.x}, y{node.y}, z{node.z}, mu{ordinate.x()};
  tep.set_variables_and_functions(
      {{"x", &x}, {"y", &y}, {"z", &z}, {"mu", &mu}});
  auto result = tep.evaluate(expression_);
  if (tep.success())
    return result;
  else
    throw std::runtime_error("Parse error at " +
                             std::to_string(tep.get_last_error_position()));
}
}  // namespace hummingbird
