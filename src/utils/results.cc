#include "utils/results.h"

#include <fstream>
#include <iomanip>
#include <limits>
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
  std::ofstream out(file_name);
  if (!out.is_open())
    throw std::runtime_error("Could not open file for writing: " + file_name);

  const size_t n_ordinates = angular_quadrature_.n_points();

  out << "node_id,x,y,z,scalar_flux";
  for (size_t i = 0; i < n_ordinates; i++)
    out << ",angular_flux_" << i << ",azimuth_" << i << ",polar_" << i;
  out << "\n";

  out << std::setprecision(std::numeric_limits<double>::max_digits10);
  for (const auto& node : nodes_) {
    out << node.id << "," << node.x << "," << node.y << "," << node.z << ","
        << node.scalar_flux;
    for (size_t i = 0; i < n_ordinates; i++) {
      const Ordinate ordinate = angular_quadrature_.GetAbscissa(i);
      out << "," << node.angular_fluxes.at(i) << "," << ordinate.azimuth()
          << "," << ordinate.polar();
    }
    out << "\n";
  }
}
}  // namespace hummingbird
