#ifndef HUMMINGBIRD_ORDINATE_H_
#define HUMMINGBIRD_ORDINATE_H_

namespace hummingbird {

/**
 * @brief Defines a solid angle, or *ordinate*, in spherical geometry
 *
 */
class Ordinate {
 public:
  /**
   * @brief Construct a new Ordinate object
   *
   * @param azimuth Azimuthal angle \f$\theta\in [-\pi,\pi]\f$
   * @param polar Polar angle \f$\varphi\in [0,\pi]\f$
   */
  Ordinate(const double azimuth, const double polar);

  /**
   * @brief Evaluate the x-direction cosine
   *
   * @return double
   */
  double XCosine();

  /**
   * @brief Evaluate the y-direction cosine
   *
   * @return double
   */
  double YCosine();

  /**
   * @brief Evaluate the z-direction cosine
   *
   * @return double
   */
  double ZCosine();

 private:
  /// @brief Azimuthal angle in radians
  const double azimuth_;

  /// @brief Polar angle in radians
  const double polar_;

  /**
   * @brief Check that input angles are within expected bounds
   *
   * @param azimuth Azimuthal angle in radians
   * @param polar Polar angle in radians
   */
  void CheckInput(const double azimuth, const double polar);
};
}  // namespace hummingbird

#endif  // HUMMINGBIRD_ORDINATE_H_
