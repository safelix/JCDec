#ifndef LOWLEVEL_HPP
#define LOWLEVEL_HPP

#include <Eigen/Dense>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/number.hpp>
#include <cmath>
#include <numeric>

// https://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
#define GCC_VERSION \
  (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

/* Type Definitions */
namespace Eigen {

template <typename ScalarType>
using VectorX = typename Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

template <typename ScalarType>
using MatrixX =
    typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

}  // namespace Eigen

namespace mp = boost::multiprecision;

/* Function Declarations and Implementations */

// MODULO:
using mp::fmod;
using std::fmod;

// https://en.cppreference.com/w/cpp/types/enable_if
template <typename Integer,
          std::enable_if_t<std::numeric_limits<Integer>::is_integer, int> = 0>
inline Integer mod(const Integer a, const Integer b) {
  return a % b;
}

// https://en.cppreference.com/w/cpp/types/enable_if
template <typename Floating,
          std::enable_if_t<!std::numeric_limits<Floating>::is_integer, int> = 0>
inline Floating mod(const Floating a, const Floating b) {
  return fmod(a, b);
}

// Round
using mp::round;
using std::round;

// ABSOLUTE VALUE:
using mp::abs;
using std::abs;

/* Fallback implementation */
template <typename Scalar>
Scalar abs(const Scalar A) {
  if (A >= 0)
    return A;
  else
    return -A;
}

// GREATEST COMMON DIVISOR
#if GCC_VERSION >= 70100
using std::gcd;
#endif
using mp::gcd;

/* Assume, only Integer Stored in Rational */
inline mp::cpp_rational gcd(const mp::cpp_rational u,
                            const mp::cpp_rational v) {
  assert(denominator(u) == 1 && denominator(v) == 1 &&
         "gcd of rationals is not defined");
  return mp::gcd(numerator(u), numerator(v));
}

/* Fallback implementation */
template <typename Scalar>
inline Scalar gcd(const Scalar U, const Scalar V) {
  Scalar u = abs(U), v = abs(V), r = 0;
  while (v != 0) {
    r = mod(u, v);
    u = v;
    v = r;
  }

  return u;
}

// INTEGER POWER:
using boost::multiprecision::pow;

/* Fallback implementation for scalar base */
template <typename ScalarB, typename ScalarE>
inline ScalarB pow(const ScalarB B, const ScalarE E) {
  if (E < 0) return 0;

  if (E == 0) return 1;

  ScalarE i = 1;
  ScalarB res = B;
  for (; i < E; i *= 2) res *= res;

  for (; i < E; i++) res *= B;

  return res;
}

/* Fallback implementation for matrix base */
template <typename ScalarB, typename ScalarE>
inline Eigen::MatrixX<ScalarB> pow(const Eigen::MatrixX<ScalarB> B,
                                   const ScalarE E) {
  eigen_assert(B.cols() == B.rows());
  Eigen::Index n = B.cols();

  if (E < 0) return Eigen::MatrixX<ScalarB>::Zero(n, n);

  if (E == 0) return Eigen::MatrixX<ScalarB>::Identity(n, n);

  ScalarE i = 1;
  Eigen::MatrixX<ScalarB> res = B;
  for (; i < E; i *= 2) res *= res;

  for (; i < E; i++) res *= B;

  return res;
}

// Numerator and denominator for rational matrices
using mp::numerator;
using mp::denominator;

inline Eigen::MatrixX<mp::cpp_int> numerator(const Eigen::MatrixXq& A) {
  return A.unaryExpr(
      [](const mp::cpp_rational frac) -> mp::cpp_int { return numerator(frac); }).cast<mp::cpp_int>(); // MSVC requires cast, gcc & clang dont
}

inline Eigen::MatrixX<mp::cpp_int> denominator(const Eigen::MatrixXq& A) {
  return A.unaryExpr(
      [](const mp::cpp_rational frac) -> mp::cpp_int { return denominator(frac); }).cast<mp::cpp_int>(); // MSVC requires cast, gcc & clang dont
}

#endif  // LOWLEVEL_HPP
