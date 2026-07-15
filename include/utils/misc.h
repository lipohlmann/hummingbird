#ifndef HUMMINGBIRD_UTILS_MISC_H
#define HUMMINGBIRD_UTILS_MISC_H

namespace hummingbird::utils {
/**
 * @brief Tests if two double-type numbers are equal using the utils::TOLERANCE
 * value. Taken from https://github.com/rcana-project/rcana/
 *
 * @param first First number to compare
 * @param second Second number to compare
 * @return true
 * @return false
 */
bool DoubleEqual(const double first, const double second);
}  // namespace hummingbird::utils

#endif  // HUMMINGBIRD_UTILS_MISC_H
