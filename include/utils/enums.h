// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_UTILS_ENUMS_H_
#define HUMMINGBIRD_UTILS_ENUMS_H_

namespace hummingbird {
/**
 * @brief Available boundary conditions
 *
 */
enum class BC { VACUUM, REFLECTIVE, NONE };

/**
 * @brief Boundary locations
 *
 */
enum class Boundary { NORTH, SOUTH, EAST, WEST };

/**
 * @brief Available volumetric source types
 *
 */
enum class SourceType { PARSED_FUNCTION };

/**
 * @brief Available finite element formulations
 *
 */
enum class FEFormulation { CG };

/**
 * @brief Available forms of the transport equation
 *
 */
enum class TransportForm { SAAF };

enum class AngularTreatmentType { DISCRETE_ORDINATES };

/**
 * @brief Angular quadrature set
 *
 */
enum class AngularQuadSet { GLT, GL, GLC };

enum class RunMode { FIXED_SOURCE };

}  // namespace hummingbird

#endif  // HUMMINGBIRD_UTILS_ENUMS_H_