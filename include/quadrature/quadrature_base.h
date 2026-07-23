#ifndef HUMMINGBIRD_QUADRATURE_QUADRATURE_BASE_H_
#define HUMMINGBIRD_QUADRATURE_QUADRATURE_BASE_H_

#include <map>
#include <vector>

namespace hummingbird::quadrature {

/**
 * @brief Simple struct holding a function value evaluated at the supplied
 * abscissa value
 *
 * @tparam T Data type of the abscissa
 */
struct QuadraturePair {
  int abscissa_index;
  double function_value;
};

/**
 * @brief Base class for defining a quadrature set
 *
 * @tparam T Type of ordinates to be used
 */
template <typename T>
class QuadratureBase {
 public:
  /**
   * @brief Construct a new Quadrature object
   *
   * @param abscissas Vector of abscissa values for the quadrature
   */
  QuadratureBase(const std::vector<T> abscissas) : abscissas_(abscissas) {};

  /**
   * @brief Get the weight value corresponding to a given abscissa value
   *
   * @param abscissa_index Abscissa index
   * @return Weight value
   */
  double GetWeight(const unsigned int abscissa_index) const {
    return weight_map_.at(abscissa_index);
  };

  /**
   * @brief Get the abscissa corresponding to the index
   *
   * @param index Abscissa index
   * @return T Abscissa object at index
   */
  T GetAbscissa(const unsigned int index) const { return abscissas_.at(index); }

  /**
   * @brief Integrate the function using the function value and its location on
   * the quadrature grid
   *
   * @param quad_pairs QuadraturePair object
   * @return Approximated integral using the quadrature set
   */
  double Integrate(const std::vector<QuadraturePair>& quad_pairs) {
    double sum = 0;
    for (auto& pair : quad_pairs) {
      auto weight = GetWeight(pair.abscissa_index);
      sum += pair.function_value * weight;
    }
    return sum;
  };

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
  /// the indices corresponding to abscissa values
  std::map<unsigned int, double> weight_map_;

  /**
   * @brief Create the weight map
   *
   */
  void CreateWeightMap() {
    int order = abscissas_.size();
    for (auto i = 0; i < order; i++) {
      double weight = ComputeWeight(i, order);
      weight_map_.insert({i, weight});
    }
  }

  virtual double ComputeWeight(const size_t k, const size_t n) = 0;
};
}  // namespace hummingbird::quadrature

#endif  // HUMMINGBIRD_QUADRATURE_QUADRATURE_BASE_H_