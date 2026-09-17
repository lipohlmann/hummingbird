#include "utils/output.h"

#include "utils/enums.h"

namespace hummingbird {
void print_header() {
  fmt::print(
      " _                               _             _     _         _ \n"
      "| |__  _   _ _ __ ___  _ __ ___ (_)_ __   __ _| |__ (_)_ __ __| |\n"
      "| '_ \\| | | | '_ ` _ \\| '_ ` _ \\| | '_ \\ / _` | '_ \\| | '__/ _` |\n"
      "| | | | |_| | | | | | | | | | | | | | | | (_| | |_) | | | | (_| |\n"
      "|_| |_|\\__,_|_| |_| |_|_| |_| |_|_|_| |_|\\__, |_.__/|_|_|  \\__,_|\n"
      "                                        |___ /                   \n "
      "\n\n");

  fmt::print(
      "    Description | Spectral Elements for Radiation Transport \n "
      "      Copyright | 2026, Liam Pohlmann\n"
      "        License | "
      "https://github.com/lipohlmann/hummingbird/blob/main/LICENSE\n\n");
}
void print_columns() {
  fmt::print(
      "   k-eff Iter.       k-eff       k-eff Error  | Scatter Iter.     "
      "Scatter L2     Scatter Error\n"
      "  =============  =============  ============= | =============  "
      "=============  =============\n");
}
void print_input_files(const std::string input, const std::string mesh) {
  fmt::print("     Input file | {}\n", input);
  fmt::print("      Mesh file | {}\n\n", mesh);
}

void print_k_status(const double k_eff, const double error,
                    const unsigned int iter) {
  if (iter > 1)
    fmt::print("  {:^13}  {:^13.5e}  {:^13.5e} |\n", iter, k_eff, error);
  else
    fmt::print("  {:^13}  {:^13.5f}  {:^13} |\n", iter, k_eff, "-");
}

void print_scatter_status(const double scatter, const double error,
                          const unsigned int iter) {
  if (iter > 1)
    fmt::print(
        "                                              | {:^13}  "
        "{:^13.5e}   {:^13.5e}\n",
        iter, scatter, error);
  else
    fmt::print(
        "                                              | {:^13}  "
        "{:^13.5e}   {:^13}\n",
        iter, scatter, "-");
}

void print_scatter_complete(const double final_k_eff, const RunMode run_mode) {
  fmt::print(
      "  "
      "------------------------------------------------------------------------"
      "-----------------"
      "\n\n");
  fmt::print("Source iterations complete.\n\n");
}
}  // namespace hummingbird
