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

NLOHMANN_JSON_SERIALIZE_ENUM(FEFormulation, {{FEFormulation::CG, "cg"}});

/**
 * @brief Available forms of the transport equation
 *
 */
enum class TransportForm { SAAF };

NLOHMANN_JSON_SERIALIZE_ENUM(TransportForm, {{TransportForm::SAAF, "saaf"}});

/**
 * @brief Angular quadrature set
 *
 */
enum class AngularQuadSet { GLT, GL, GLC };

NLOHMANN_JSON_SERIALIZE_ENUM(AngularQuadSet,
                             {{AngularQuadSet::GLT, "gauss_legendre_trapezoid"},
                              {AngularQuadSet::GLC, "gauss_legendre_chebyshev"},
                              {AngularQuadSet::GL, "gauss_legendre"}});

/**
 * @brief Run modes available
 *
 */
enum class RunMode { FIXED_SOURCE };

NLOHMANN_JSON_SERIALIZE_ENUM(RunMode,
                             {{RunMode::FIXED_SOURCE, "fixed_source"}});

/**
 * @brief Output forms available for export
 *
 */
enum class OutputFormat { VTK, CSV };

NLOHMANN_JSON_SERIALIZE_ENUM(OutputFormat, {{OutputFormat::CSV, "csv"},
                                            {OutputFormat::VTK, "vtk"}});

}  // namespace hummingbird

#endif  // HUMMINGBIRD_UTILS_ENUMS_H_
