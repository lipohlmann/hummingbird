#ifndef HUMMINGBIRD_MATH_LEGENDRE_POLYNOMIALS_H
#define HUMMINGBIRD_MATH_LEGENDRE_POLYNOMIALS_H

namespace hummingbird::math {
/**
 * @brief Computes the Legendre polynomial of degree n at a point x using the
 * recursion formula. Legendre polynomials are defined on the interval \f$x\in
 * [-1,1]\f$.
 *
 * @param n Order of polynomial to evaluate
 * @param x Point to evaluate at
 * @return double
 */
double LegendrePolynomial(const int n, const double x);

}  // namespace hummingbird::math

#endif  // HUMMINGBIRD_MATH_LEGENDRE_POLYNOMIALS_H
