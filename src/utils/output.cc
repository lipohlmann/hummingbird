#include "output.h"

#include "enums.h"

namespace starling {
void print_header() {
  fmt::print(
      "                                       ###                              "
      "  \n"
      "                                        ###    #                        "
      "  \n"
      "               #                         ##   ###                       "
      "  \n"
      "              ##                         ##    #                        "
      "  \n"
      "              ##                         ##                             "
      "  \n"
      "   /###     ######## /###   ###  /###    ##  ###   ###  /###     /###   "
      "  \n"
      "  / #### / ######## / ###  / ###/ #### / ##   ###   ###/ #### / /  ###  "
      "/ \n"
      " ##  ###/     ##   /   ###/   ##   ###/  ##    ##    ##   ###/ /    "
      "###/  \n"
      "####          ##  ##    ##    ##         ##    ##    ##    ## ##     ## "
      "  \n"
      "  ###         ##  ##    ##    ##         ##    ##    ##    ## ##     ## "
      "  \n"
      "    ###       ##  ##    ##    ##         ##    ##    ##    ## ##     ## "
      "  \n"
      "      ###     ##  ##    ##    ##         ##    ##    ##    ## ##     ## "
      "  \n"
      " /###  ##     ##  ##    /#    ##         ##    ##    ##    ## ##     ## "
      "  \n"
      "/ #### /      ##   ####/ ##   ###        ### / ### / ###   ### ######## "
      "  \n"
      "   ###/        ##   ###   ##   ###        ##/   ##/   ###   ###  ### "
      "###  \n"
      "                                                                      "
      "### \n"
      "                                                                ####   "
      "###\n"
      "                                                              /######  "
      "/# \n"
      "                                                             /     ###/ "
      "  \n\n");

  fmt::print(
      "    Description | Radiance Cascades on Voxel Geometry for "
      "Nuclear "
      "Reactor Physics\n"
      "      Copyright | 2026, Liam Pohlmann\n"
      "        License | "
      "https://github.com/lipohlmann/starling/blob/main/LICENSE\n\n");
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

void print_scatter_complete(const double final_k_eff,
                            const SimulationType sim_type) {
  fmt::print(
      "  "
      "------------------------------------------------------------------------"
      "-----------------"
      "\n\n");
  fmt::print("Source iterations complete.\n\n");

  if (sim_type == SimulationType::K_EIGENVALUE)
    fmt::print("Final k-eff | {:.5f}\n\n", final_k_eff);
}
}  // namespace starling
