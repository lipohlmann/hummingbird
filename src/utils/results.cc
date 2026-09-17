#include "utils/results.h"

#include <stdexcept>

namespace hummingbird {
Results::Results(const std::string& name, const OutputFormat output_format,
                 const std::vector<Node>& nodes)
    : name_(name), output_format_(output_format), nodes_(nodes) {}

void Results::Export() {
  switch (output_format_) {
    case OutputFormat::CSV:
      ToCSV();
      break;

    default:
      throw std::runtime_error("Output format not yet supported.");
      break;
  }
}
}  // namespace hummingbird
