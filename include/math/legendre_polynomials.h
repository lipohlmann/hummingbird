/* The following was written heavily based on the teachings and codes found in
 * the book Numerical Methods in Physics with Python, Second Edition by Alex
 * Gezerlis. For more information on this text, please see
 * https://numphyspy.org/. */

#ifndef HUMMINGBIRD_MATH_LEGENDRE_POLYNOMIALS_H
#define HUMMINGBIRD_MATH_LEGENDRE_POLYNOMIALS_H

#include <vector>

namespace hummingbird::math {
/**
 * @brief Computes the Legendre polynomial of degree n at a point x using the
 * formula:
 *
 * \f[
 * P_{n+1}(x)=\frac{(2n+1)xP_{n}(x)-nP_{n-1}(x)}{n+1}
 * \f]
 *
 * with:
 *
 * \f[
 * P_0=1,\quad P_1=x
 * \f]
 *
 * The above is implemented to evaluate the n-th order Legendre polynomial by
 * building from the bottom up, as opposed to using a recursion relation. This
 * greatly improves the evaluation speed. This algorithm was taken from
 * https://github.com/CambridgeUniversityPress/NumericalMethodsPhysicsWithPython/blob/master/second_edition/codes/legendre.py
 * written by Alex Gezerlis for his book, Numerical Methods in Physics With
 * Python, Second Edition.
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
 * The value of the derivative at the bounds is special:
 *
 * \f[
 * P'_n(x)=\begin{cases} (-1)^{n-1}\frac{n(n+1)}{2} & x=-1
 * \\ \frac{n(n+1)}{2} &x=1 \end{cases}
 * \f]
 *
 * @param n
 * @param x
 * @return double
 */
double LegendrePolynomialPrime(const int n, const double x);

double LegendrePolynomialPrimePrime(const int n, const double x);

/**
 * @brief Compute all roots of the Legendre polynomial of order n
 *
 * @param n Order
 * @return std::vector<double>
 */
std::vector<double> AllLegendreRoots(const int n);

/**
 * @brief Computes the k-th root of the n-th order Legendre polynomial by first
 * approximating the root with
 *
 * \f[
 * x_k^{(0)}\approx \cos \left(\frac{4k+3}{4n+2}\pi\right)
 * \f]
 * then using Newton's method to converge to the root.
 *
 * @param n Order of the polynomial
 * @param k Root to find
 * @return double
 */
double LegendreRoot(const int n, const int k);

double LegendrePrimeRoot(const int n, const int k);

}  // namespace hummingbird::math

#endif  // HUMMINGBIRD_MATH_LEGENDRE_POLYNOMIALS_H
