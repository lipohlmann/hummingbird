#include "quadrature/quadrature.h"

namespace hummingbird::quadrature {
template <typename T>
Quadrature<T>::Quadrature(const std::vector<T> abscissas)
    : abscissas_(abscissas) {}

template <typename T>
double Quadrature<T>::Integrate(
    const std::vector<QuadraturePair<T>>& quad_pairs) {
  double sum = 0;
  for (auto& pair : quad_pairs) {
    auto weight = GetWeight(pair.abscissa);
    sum += pair.function_value * weight;
  }
  return sum;
}
}  // namespace hummingbird::quadrature
