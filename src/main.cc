#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

#include "input_parameters.h"
#include "utils/json.h"

using nlohmann::json;
using namespace hummingbird;

int main(int argc, char **argv) {
  if (argc < 2) {
    printf("Usage: %s <input_file_name>\n", argv[0]);
    return 1;
  };

  std::filesystem::path working_directory(argv[1]);
  working_directory.remove_filename();

  json user_input_json = JSONFromFile(argv[1]);

  InputParams input_params = user_input_json.get<InputParams>();

  // build Mesh object

  // form local matrices

  // assemble to global system

  // solve system (should be a one-liner?)

  // check source iteration convergence

  // export results

  return 0;
}
