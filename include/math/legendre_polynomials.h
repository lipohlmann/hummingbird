#ifndef HUMMINGBIRD_MATH_LEGENDRE_POLYNOMIALS_H
#define HUMMINGBIRD_MATH_LEGENDRE_POLYNOMIALS_H

namespace hummingbird::math {
/**
 * @brief Computes the Legendre polynomial of degree n at a point x using the
 * recursion formula:
 *
 * \f[
 * P_{n+1}(x)=\frac{(2n+1)xP_{n}(x)-nP_{n-1}(x)}{n+1}
 * \f]
 *
 *  Legendre polynomials are defined on the interval \f$x\in
 * [-1,1]\f$.
 *
 * @param n Order of polynomial to evaluate
 * @param x Point to evaluate at
 * @return double
 */
double LegendrePolynomial(const int n, const double x);

/**
 * @brief Computes the first derivative of the Legendre polynomial of degree n
 * at a point x using the recurrence relation:
 *
 * \f[
 * P'_n (x)=\frac{nP_{n-1}(x)-nxP_n(x)}{1-x^2}
 * \f]
 *
 * @param n
 * @param x
 * @return double
 */
double LegendrePolynomialDerivative(const int n, const double x);

}  // namespace hummingbird::math

#endif  // HUMMINGBIRD_MATH_LEGENDRE_POLYNOMIALS_H
