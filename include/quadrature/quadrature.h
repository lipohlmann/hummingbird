#ifndef HUMMINGBIRD_QUADRATURE_QUADRATURE_H_
#define HUMMINGBIRD_QUADRATURE_QUADRATURE_H_

#include <map>
#include <vector>

namespace hummingbird::quadrature {

/**
 * @brief Simple struct holding a function value evaluated at the supplied
 * abscissa value
 *
 * @tparam T Data type of the abscissa
 */
template <typename T>
struct QuadraturePair {
  T abscissa;
  double function_value;
};

/**
 * @brief Base class for defining a quadrature set
 *
 * @tparam T Type of ordinates to be used
 */
template <typename T>
class Quadrature {
 public:
  Quadrature(const std::vector<T> abscissas);

  /**
   * @brief Get the weight value corresponding to a given abscissa value
   *
   * @param abscissa Abscissa value
   * @return Weight value
   */
  double GetWeight(const T& abscissa) const {
    return weight_map_.get(abscissa);
  };

  /**
   * @brief Integrate the function using the function value and its location on
   * the quadrature grid
   *
   * @param quad_pairs QuadraturePair object
   * @return Approximated integral using the quadrature set
   */
  double Integrate(const std::vector<QuadraturePair<T>>& quad_pairs);

  /**
   * @brief Getter for the abscissas in the quadrature set
   *
   * @return const std::vector<T>&
   */
  const std::vector<T>& abscissas() const { return abscissas_; }

 protected:
  /// @brief Abscissa values
  const std::vector<T> abscissas_;

  /// @brief Map to store the weights of the quadrature set where the keys are
  /// their corresponding ordinates
  std::map<T, double> weight_map_;
};
}  // namespace hummingbird::quadrature

#endif  // HUMMINGBIRD_QUADRATURE_QUADRATURE_H_