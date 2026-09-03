// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Liam Pohlmann

#ifndef HUMMINGBIRD_ENUMS_H_
#define HUMMINGBIRD_ENUMS_H_

namespace hummingbird {
enum class BC { TRANSMISSIVE, REFLECTIVE, NONE };
enum class Boundary { NORTH, SOUTH, EAST, WEST };
enum class SourceType { MANUFACTURED, CONSTANT };
enum class Formulation { CG };

enum class AngularQuadSet { GLT, GL, GLC };

}  // namespace hummingbird

#endif  // HUMMINGBIRD_ENUMS_H_