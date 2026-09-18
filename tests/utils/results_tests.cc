#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "mesh/node.h"
#include "quadrature/angular/gauss_legendre_trapezoid.h"
#include "quadrature/angular/ordinate.h"
#include "utils/enums.h"
#include "utils/results.h"

namespace hummingbird {

namespace {

// Helper: reads back the `<name>_results.csv` file produced by ToCSV, and
// removes it once the test is done so runs don't leave stray output behind.
class ResultsFile {
 public:
  explicit ResultsFile(const std::string& name)
      : path_(name + "_results.csv") {}

  ~ResultsFile() { std::filesystem::remove(path_); }

  std::vector<std::string> ReadLines() const {
    std::ifstream in(path_);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) lines.push_back(line);
    return lines;
  }

 private:
  std::filesystem::path path_;
};

// Formats a value the same way Results::ToCSV does, so expected rows can be
// built without duplicating the exact precision/format used by ToCSV.
template <typename T>
std::string FormatValue(const T& value) {
  std::ostringstream out;
  out << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
  return out.str();
}

}  // namespace

TEST(ResultsToCSVTest, WritesHeaderAndRowForSingleOrdinateSingleNode) {
  GaussLegendreTrapezoid quad(1, 1);

  Node node(0, 1.0, 2.0, 3.0);
  node.scalar_flux = 5.0;
  node.angular_fluxes = {7.5};
  std::vector<Node> nodes = {node};

  const std::string name = "results_test_single";
  ResultsFile file(name);
  Results results(name, OutputFormat::CSV, nodes, quad, /*dimension=*/1);
  results.Export();

  std::vector<std::string> lines = file.ReadLines();
  ASSERT_EQ(lines.size(), 2u);
  EXPECT_EQ(lines[0], "node_id,x,y,z,scalar_flux,angular_flux_0,"
                      "direction_cosine_0");

  const Ordinate ordinate = quad.GetAbscissa(0);
  const std::string expected_row =
      FormatValue(node.id) + "," + FormatValue(node.x) + "," +
      FormatValue(node.y) + "," + FormatValue(node.z) + "," +
      FormatValue(node.scalar_flux) + "," +
      FormatValue(node.angular_fluxes[0]) + "," + FormatValue(ordinate.x());
  EXPECT_EQ(lines[1], expected_row);
}

TEST(ResultsToCSVTest, WritesColumnGroupPerOrdinateAndRowPerNode) {
  GaussLegendreTrapezoid quad(2, 2);
  const size_t n_ordinates = quad.n_points();
  ASSERT_GT(n_ordinates, 1u);

  Node node_0(0, 0.0, 0.0, 0.0);
  node_0.scalar_flux = 1.0;
  node_0.angular_fluxes.resize(n_ordinates);
  for (size_t i = 0; i < n_ordinates; i++)
    node_0.angular_fluxes[i] = static_cast<double>(i) + 1.0;

  Node node_1(1, 4.0, 5.0, 6.0);
  node_1.scalar_flux = 2.0;
  node_1.angular_fluxes.resize(n_ordinates);
  for (size_t i = 0; i < n_ordinates; i++)
    node_1.angular_fluxes[i] = 2.0 * (static_cast<double>(i) + 1.0);

  std::vector<Node> nodes = {node_0, node_1};

  const std::string name = "results_test_multi";
  ResultsFile file(name);
  Results results(name, OutputFormat::CSV, nodes, quad, /*dimension=*/1);
  results.Export();

  std::vector<std::string> lines = file.ReadLines();
  ASSERT_EQ(lines.size(), 3u);  // header + one row per node

  std::string expected_header = "node_id,x,y,z,scalar_flux";
  for (size_t i = 0; i < n_ordinates; i++)
    expected_header += ",angular_flux_" + std::to_string(i) +
                       ",direction_cosine_" + std::to_string(i);
  EXPECT_EQ(lines[0], expected_header);

  const std::vector<Node> expected_nodes = {node_0, node_1};
  for (size_t row = 0; row < expected_nodes.size(); row++) {
    const Node& node = expected_nodes[row];
    std::string expected_row = FormatValue(node.id) + "," +
                               FormatValue(node.x) + "," + FormatValue(node.y) +
                               "," + FormatValue(node.z) + "," +
                               FormatValue(node.scalar_flux);
    for (size_t i = 0; i < n_ordinates; i++) {
      const Ordinate ordinate = quad.GetAbscissa(i);
      expected_row += "," + FormatValue(node.angular_fluxes[i]) + "," +
                      FormatValue(ordinate.x());
    }
    EXPECT_EQ(lines[row + 1], expected_row);
  }
}

TEST(ResultsToCSVTest, NonOneDimensionalUsesAzimuthAndPolarColumns) {
  GaussLegendreTrapezoid quad(1, 1);

  Node node(0, 1.0, 2.0, 3.0);
  node.scalar_flux = 5.0;
  node.angular_fluxes = {7.5};
  std::vector<Node> nodes = {node};

  const std::string name = "results_test_non_1d";
  ResultsFile file(name);
  Results results(name, OutputFormat::CSV, nodes, quad, /*dimension=*/2);
  results.Export();

  std::vector<std::string> lines = file.ReadLines();
  ASSERT_EQ(lines.size(), 2u);
  EXPECT_EQ(lines[0],
            "node_id,x,y,z,scalar_flux,angular_flux_0,azimuth_0,polar_0");

  const Ordinate ordinate = quad.GetAbscissa(0);
  const std::string expected_row =
      FormatValue(node.id) + "," + FormatValue(node.x) + "," +
      FormatValue(node.y) + "," + FormatValue(node.z) + "," +
      FormatValue(node.scalar_flux) + "," +
      FormatValue(node.angular_fluxes[0]) + "," +
      FormatValue(ordinate.azimuth()) + "," + FormatValue(ordinate.polar());
  EXPECT_EQ(lines[1], expected_row);
}

TEST(ResultsToCSVTest, ThrowsWhenAngularFluxesShorterThanQuadrature) {
  GaussLegendreTrapezoid quad(2, 2);

  Node node(0, 0.0, 0.0, 0.0);
  node.scalar_flux = 0.0;
  node.angular_fluxes = {1.0};  // fewer entries than quad.n_points()
  std::vector<Node> nodes = {node};

  const std::string name = "results_test_short_fluxes";
  ResultsFile file(name);
  Results results(name, OutputFormat::CSV, nodes, quad, /*dimension=*/1);

  EXPECT_THROW(results.Export(), std::out_of_range);
}

}  // namespace hummingbird
