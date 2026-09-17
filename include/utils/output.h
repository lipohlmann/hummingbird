// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_OUTPUT_H_
#define HUMMINGBIRD_OUTPUT_H_

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <string>

#include "utils/enums.h"

namespace hummingbird {

void print_header();
void print_input_files(const std::string input, const std::string mesh);
void print_columns();
void print_k_status(const double k_eff, const double error,
                    const unsigned int iter);
void print_scatter_status(const double scatter, const double error,
                          const unsigned int iter);
void print_scatter_complete(const double final_k_eff, const RunMode sim_typ);

}  // namespace hummingbird

#endif  // HUMMINGBIRD_OUTPUT_H_
