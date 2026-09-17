// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_UTILS_RESULTS_H_
#define HUMMINGBIRD_UTILS_RESULTS_H_

#include <string>
#include <vector>

#include "mesh/node.h"
#include "quadrature/angular/ordinate.h"
#include "quadrature/quadrature_base.h"
#include "utils/enums.h"

namespace hummingbird {

/**
 * @brief Handles exporting simulation results to a file
 *
 */
class Results {
 public:
  /**
   * @brief Construct a new Results object
   *
   * @param name Problem name, used as the output file name prefix
   * @param output_format Desired output file format
   * @param nodes Mesh nodes containing the solution to export
   * @param angular_quadrature Angular quadrature set used to solve the
   * problem. Needed to report the azimuth and polar angles corresponding to
   * each angular flux value
   */
  Results(const std::string& name, const OutputFormat output_format,
          const std::vector<Node>& nodes,
          const QuadratureBase<Ordinate>& angular_quadrature);

  /**
   * @brief Export the results to a file in the requested output format
   *
   */
  void Export();

 private:
  /// @brief Problem name, used as the output file name prefix
  const std::string& name_;

  /// @brief Desired output file format
  const OutputFormat output_format_;

  /// @brief Mesh nodes containing the solution to export
  const std::vector<Node>& nodes_;

  /// @brief Angular quadrature set used to solve the problem
  const QuadratureBase<Ordinate>& angular_quadrature_;

  /**
   * @brief Write the results to a CSV file named "<name>_results.csv". Each
   * row holds a node's ID, coordinates, and scalar flux, followed by the
   * azimuth angle, polar angle, and angular flux value for each ordinate in
   * the angular quadrature set
   *
   */
  void ToCSV();
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_UTILS_RESULTS_H_
