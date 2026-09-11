#include "utils/json.h"

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <stdexcept>

using nlohmann::json;

namespace hummingbird {
json JSONFromFile(const std::string filename) {
  std::filesystem::path file_path = filename;
  std::ifstream fs_file(filename);
  if (!fs_file.is_open()) {
    std::filesystem::path absolute_path = std::filesystem::absolute(file_path);
    throw std::runtime_error("Could not open mesh file: " +
                             absolute_path.string());
  }

  return json::parse(fs_file);
}
}  // namespace hummingbird
