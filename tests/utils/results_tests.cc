#include "utils/results.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "mesh/node.h"
#include "quadrature/angular/gauss_legendre_trapezoid.h"
#include "quadrature/angular/ordinate.h"
#include "utils/enums.h"

namespace hummingbird {

namespace {

// Reads back every line of a file into a vector of strings.
std::vector<std::string> ReadLines(const std::string& path) {
  std::ifstream file(path);
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(file, line)) lines.push_back(line);
  return lines;
}

}  // namespace

TEST(ResultsToCSVTest, WritesHeaderAndSingleOrdinateRow) {
  const std::string name = "results_test_single_ordinate";
  const std::string path = name + "_results.csv";
  std::filesystem::remove(path);

  GaussLegendreTrapezoid quad(1, 1);
  const Ordinate ordinate = quad.GetAbscissa(0);

  Node node(0, 1.0, 2.0, 3.0);
  node.scalar_flux = 4.0;
  node.angular_fluxes = {5.0};

  std::vector<Node> nodes = {node};
  Results results(name, OutputFormat::CSV, nodes, quad);
  results.Export();

  ASSERT_TRUE(std::filesystem::exists(path));
  std::vector<std::string> lines = ReadLines(path);
  ASSERT_EQ(lines.size(), 2u);
  EXPECT_EQ(lines[0],
           "node_id,x,y,z,scalar_flux,angular_flux_0,azimuth_0,polar_0");

  std::ostringstream expected_row;
  expected_row << std::setprecision(15) << "0,1,2,3,4,5," << ordinate.azimuth()
              << "," << ordinate.polar();
  EXPECT_EQ(lines[1], expected_row.str());

  std::filesystem::remove(path);
}

TEST(ResultsToCSVTest, WritesOneColumnGroupPerOrdinateForMultipleNodes) {
  const std::string name = "results_test_multiple_ordinates";
  const std::string path = name + "_results.csv";
  std::filesystem::remove(path);

  GaussLegendreTrapezoid quad(2, 2);

  Node node_1(0, 0.0, 0.0, 0.0);
  node_1.scalar_flux = 1.0;
  node_1.angular_fluxes.resize(quad.n_points());
  for (size_t i = 0; i < quad.n_points(); i++)
    node_1.angular_fluxes[i] = static_cast<double>(i);

  Node node_2(1, 1.0, 1.0, 1.0);
  node_2.scalar_flux = 2.0;
  node_2.angular_fluxes.resize(quad.n_points());
  for (size_t i = 0; i < quad.n_points(); i++)
    node_2.angular_fluxes[i] = static_cast<double>(i) + 10.0;

  std::vector<Node> nodes = {node_1, node_2};
  Results results(name, OutputFormat::CSV, nodes, quad);
  results.Export();

  ASSERT_TRUE(std::filesystem::exists(path));
  std::vector<std::string> lines = ReadLines(path);
  ASSERT_EQ(lines.size(), 3u);

  std::ostringstream expected_header;
  expected_header << "node_id,x,y,z,scalar_flux";
  for (size_t i = 0; i < quad.n_points(); i++)
    expected_header << ",angular_flux_" << i << ",azimuth_" << i << ",polar_"
                    << i;
  EXPECT_EQ(lines[0], expected_header.str());

  for (const auto& node : {node_1, node_2}) {
    std::ostringstream expected_row;
    expected_row << std::setprecision(15) << node.id << "," << node.x << ","
                << node.y << "," << node.z << "," << node.scalar_flux;
    for (size_t i = 0; i < quad.n_points(); i++) {
      const Ordinate ordinate = quad.GetAbscissa(i);
      expected_row << "," << node.angular_fluxes.at(i) << ","
                  << ordinate.azimuth() << "," << ordinate.polar();
    }
    EXPECT_EQ(lines[node.id + 1], expected_row.str());
  }

  std::filesystem::remove(path);
}

}  // namespace hummingbird
