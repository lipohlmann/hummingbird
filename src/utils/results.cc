#include "utils/results.h"

#include <fstream>
#include <stdexcept>

namespace hummingbird {
Results::Results(const std::string& name, const OutputFormat output_format,
                 const std::vector<Node>& nodes,
                 const QuadratureBase<Ordinate>& angular_quadrature)
    : name_(name),
      output_format_(output_format),
      nodes_(nodes),
      angular_quadrature_(angular_quadrature) {}

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

void Results::ToCSV() {
  const std::string file_name = name_ + "_results.csv";
  std::ofstream file(file_name);
  if (!file.is_open())
    throw std::runtime_error("Unable to open file for writing: " + file_name);

  const size_t n_ordinates = angular_quadrature_.n_points();

  file << "node_id,x,y,z,scalar_flux";
  for (size_t n = 0; n < n_ordinates; n++)
    file << ",azimuth_" << n << ",polar_" << n << ",angular_flux_" << n;
  file << "\n";

  for (const auto& node : nodes_) {
    file << node.id << "," << node.x << "," << node.y << "," << node.z << ","
        << node.scalar_flux;
    for (size_t n = 0; n < n_ordinates; n++) {
      const auto ordinate = angular_quadrature_.GetAbscissa(n);
      file << "," << ordinate.azimuth() << "," << ordinate.polar() << ","
          << node.angular_fluxes.at(n);
    }
    file << "\n";
  }
}
}  // namespace hummingbird
