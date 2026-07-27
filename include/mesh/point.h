#ifndef HUMMINGBIRD_MESH_POINT_H_
#define HUMMINGBIRD_MESH_POINT_H_

namespace hummingbird::mesh {
/**
 * @brief Defines a point in Cartesian space
 *
 */
struct Point {
  double x = 0;
  double y = 0;
  double z = 0;
};

}  // namespace hummingbird::mesh

#endif  // HUMMINGBIRD_MESH_POINT_H_