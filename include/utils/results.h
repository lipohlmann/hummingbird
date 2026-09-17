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
class Results {
 public:
  Results(const std::string& name, const OutputFormat output_format,
          const std::vector<Node>& nodes,
          const QuadratureBase<Ordinate>& angular_quadrature);

  void Export();

 private:
  const std::string& name_;
  const OutputFormat output_format_;
  const std::vector<Node>& nodes_;
  const QuadratureBase<Ordinate>& angular_quadrature_;

  void ToCSV();
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_UTILS_RESULTS_H_
