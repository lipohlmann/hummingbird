// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef STARLING_UTILS_JSON_H_
#define STARLING_UTILS_JSON_H_

#include <nlohmann/json.hpp>
#include <string>

using nlohmann::json;

namespace starling::utils {

/**
 * @brief Create a json object from a string. Will throw an error if the file
 * could not be opened
 *
 * @param filename File name (relative path) to be opened
 * @return json
 */
json JSONFromFile(const std::string filename);

}  // namespace starling::utils

#endif  // STARLING_UTILS_JSON_H_
