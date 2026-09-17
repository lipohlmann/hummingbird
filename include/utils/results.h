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
 * @brief Exports simulation results to a file
 *
 */
class Results {
 public:
  /**
   * @brief Construct a new Results object
   *
   * @param name Problem name, used to build the output file name
   * @param output_format Format to export results in
   * @param nodes Mesh nodes holding the solution values to export
   */
  Results(const std::string& name, const OutputFormat output_format,
          const std::vector<Node>& nodes,
          const QuadratureBase<Ordinate>& angular_quadrature);

  /**
   * @brief Export the results in the format given by output_format_
   *
   * @throw std::runtime_error if output_format_ is not yet supported
   */
  void Export();

 private:
  /// @brief Problem name, used to build the output file name
  const std::string& name_;

  /// @brief Format to export results in
  const OutputFormat output_format_;

  /// @brief Mesh nodes holding the solution values to export
  const std::vector<Node>& nodes_;

  /// @brief Angular quadrature set
  const QuadratureBase<Ordinate>& angular_quadrature_;

  /**
   * @brief Export the results as a CSV file
   *
   */
  void ToCSV();
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_UTILS_RESULTS_H_
