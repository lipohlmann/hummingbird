#include "utils/results.h"

#include <fstream>
#include <iomanip>
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
  std::ofstream file(name_ + "_results.csv");
  if (!file.is_open())
    throw std::runtime_error("Unable to open file for CSV export: " + name_ +
                             "_results.csv");
  file << std::setprecision(15);

  file << "node_id,x,y,z,scalar_flux";
  for (size_t n = 0; n < angular_quadrature_.n_points(); n++)
    file << ",angular_flux_" << n << ",azimuth_" << n << ",polar_" << n;
  file << "\n";

  for (const auto& node : nodes_) {
    file << node.id << "," << node.x << "," << node.y << "," << node.z << ","
        << node.scalar_flux;
    for (size_t n = 0; n < angular_quadrature_.n_points(); n++) {
      const auto& ordinate = angular_quadrature_.GetAbscissa(n);
      file << "," << node.angular_fluxes.at(n) << "," << ordinate.azimuth()
          << "," << ordinate.polar();
    }
    file << "\n";
  }
}
}  // namespace hummingbird
