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

/**
 * @brief Computes the second derivative of the Legendre polynomial of degree n
 * using Legendre's differential equation:
 *
 * \f[
 * (1-x^2)P_n''(x)-2xP_n'(x)+n(n+1)P_n(x)=0
 * \f]
 *
 * @param n Order
 * @param x Location to evaluate at
 * @return double
 */
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
 * approximating the root with ApproximateLegendreRoot, then using Newton's
 * method to converge to the root.
 *
 * @param n Order of the polynomial
 * @param k Root to find
 * @return double
 */
double LegendreRoot(const int n, const int k);

/**
 * @brief Approximates the k-th root of the n-th order Legendre polynomial with:
 *
 * \f[
 * x_k^{(0)}\approx \cos \left(\frac{4k+3}{4n+2}\pi\right)
 * \f]
 *
 * @param n Order
 * @param k Root number
 * @return double
 * @throw std::invalid_argument Root number must be less than polynomial order
 */
double ApproximateLegendreRoot(const int n, const int k);

/**
 * @brief Computes all roots of the \f$P'_n(x)\f$ polynomial. A total of n-1
 * roots will be computed and returned
 *
 * @param n Order
 * @return std::vector<double>
 */
std::vector<double> AllLegendrePrimeRoots(const int n);

/**
 * @brief Computes the k-th root of the n-th order first derivative of the
 * Legendre polynomial by first approximating with ApproximateLegendrePrimeRoot,
 * then using Newton's method to converge
 *
 * @param n Order
 * @param k Root number
 * @return double
 */
double LegendrePrimeRoot(const int n, const int k);

/**
 * @brief Approximates the k-th root of the first derivative of the n-th degree
 * Legendre polynomial using the average of the approximations
 * (ApproximateLegendreRoot()) of the surrounding roots of the n-th degree
 * Legendre polynomial. Note that zero-indexing is used, so the roots of the
 * \f$P'_4 (x)\f$ polynomial are \f$k=0,1,2\f$. The approximations of these are
 * given by:
 *
 * \f[
 * x_k\approx \frac{1}{2} \left( \cos \left(\frac{4k+3}{4n+2}\pi\right) + \cos
 * \left(\frac{4(k+1)+3}{4n+2}\pi\right) \right) \f]
 *
 * @param n Order
 * @param k Root number (zero-indexed)
 * @return double
 */
double ApproximateLegendrePrimeRoot(const int n, const int k);

}  // namespace hummingbird::math

#endif  // HUMMINGBIRD_MATH_LEGENDRE_POLYNOMIALS_H
