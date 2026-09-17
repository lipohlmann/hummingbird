// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef STARLING_OUTPUT_H_
#define STARLING_OUTPUT_H_

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <string>

#include "enums.h"

namespace starling {

void print_header();
void print_input_files(const std::string input, const std::string mesh);
void print_columns();
void print_k_status(const double k_eff, const double error,
                    const unsigned int iter);
void print_scatter_status(const double scatter, const double error,
                          const unsigned int iter);
void print_scatter_complete(const double final_k_eff,
                            const SimulationType sim_typ);

}  // namespace starling

#endif  // STARLING_OUTPUT_H_
