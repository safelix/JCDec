#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <Eigen/Dense>

/* Type Definitions */
namespace Eigen {

template <typename ScalarType>
using VectorX = typename Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

template <typename ScalarType>
using MatrixX =
    typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

}  // namespace Eigen

template <typename Scalar>
using enable_if_integer =
    typename std::enable_if_t<std::numeric_limits<Scalar>::is_integer, int>;

template <typename Scalar>
using enable_if_not_integer =
    typename std::enable_if_t<!std::numeric_limits<Scalar>::is_integer, int>;

/* Function Definitions */

template <typename Scalar>
inline bool poly_isConst(const Eigen::VectorX<Scalar> &P);

template <typename Scalar>
inline bool poly_isZero(const Eigen::VectorX<Scalar> &P);

template <typename Scalar>
inline bool poly_isOne(const Eigen::VectorX<Scalar> &P);

template <typename Scalar>
inline bool poly_isEqual(const Eigen::VectorX<Scalar> &P,
                         const Eigen::VectorX<Scalar> &Q);

template <typename Scalar>
inline Eigen::VectorX<Scalar> poly_add(const Eigen::VectorX<Scalar> &P,
                                       const Eigen::VectorX<Scalar> &Q);

template <typename Scalar>
inline Eigen::VectorX<Scalar> poly_sub(const Eigen::VectorX<Scalar> &P,
                                       const Eigen::VectorX<Scalar> &Q);

template <typename Scalar>
inline Eigen::VectorX<Scalar> poly_mult(const Eigen::VectorX<Scalar> &P,
                                        const Eigen::VectorX<Scalar> &Q);

/* Polynomial pesudo division for integers: c*U(X) = V(X)q(X) + r(X)
    => 'bool pseudo' has no effect, but is needed for overloading */
template <typename Scalar, enable_if_integer<Scalar> = 0>
inline Eigen::VectorX<Scalar> poly_div(const Eigen::VectorX<Scalar> &P,
                                       const Eigen::VectorX<Scalar> &Q,
                                       const bool pseudo = true);

/* Polynomial (pseudo-) division for arbitrary types */
template <typename Scalar, enable_if_not_integer<Scalar> = 0>
inline Eigen::VectorX<Scalar> poly_div(const Eigen::VectorX<Scalar> &P,
                                       const Eigen::VectorX<Scalar> &Q,
                                       const bool pseudo = false);

/* Polynomial pseudo modulo for integers: c*U(X) = V(X)q(X) + r(X)
    => 'bool pseudo' has no effect, but is needed for overloading */
template <typename Scalar, enable_if_integer<Scalar> = 0>
inline Eigen::VectorX<Scalar> poly_mod(const Eigen::VectorX<Scalar> &P,
                                       const Eigen::VectorX<Scalar> &Q,
                                       const bool pseudo = true);

/* Polynomial (pseudo-) modulo for arbitrary types */
template <typename Scalar, enable_if_not_integer<Scalar> = 0>
inline Eigen::VectorX<Scalar> poly_mod(const Eigen::VectorX<Scalar> &P,
                                       const Eigen::VectorX<Scalar> &Q,
                                       const bool pseudo = false);

/* Polynomial pseudo quotient and remainder c*U(X) = V(X)q(X) + r(X)
    => returns q(x), stores r(x) in the VectorX<Scalar> &rem
    => 'bool pseudo' has no effect, but is needed for overloading */
template <typename Scalar, enable_if_integer<Scalar> = 0>
inline Eigen::VectorX<Scalar> poly_quorem(const Eigen::VectorX<Scalar> &P,
                                          const Eigen::VectorX<Scalar> &Q,
                                          Eigen::VectorX<Scalar> &rem,
                                          const bool pseudo = true);

/* Polynomial (pseudo-) quotient and remainder
    => returns q(x), stores r(x) in the VectorX<Scalar> &rem */
template <typename Scalar, enable_if_not_integer<Scalar> = 0>
inline Eigen::VectorX<Scalar> poly_quorem(const Eigen::VectorX<Scalar> &P,
                                          const Eigen::VectorX<Scalar> &Q,
                                          Eigen::VectorX<Scalar> &rem,
                                          const bool pseudo = false);

/* Polynomial greatest commmon divisor, using the resultant algorithm
    => 'bool resultant' has no effect, but is needed for overloading */
template <typename Scalar, enable_if_integer<Scalar> = 0>
inline Eigen::VectorX<Scalar> poly_gcd(const Eigen::VectorX<Scalar> &P,
                                       const Eigen::VectorX<Scalar> &Q,
                                       const bool resultant = true);

/* Polynomial greatest commmon divisor, using the gcd or the resultant algorithm
 */
template <typename Scalar, enable_if_not_integer<Scalar> = 0>
inline Eigen::VectorX<Scalar> poly_gcd(const Eigen::VectorX<Scalar> &P,
                                       const Eigen::VectorX<Scalar> &Q,
                                       const bool resultant = false);

/* Polynomial modular inverse, coefficients could grow exponentially */
template <typename Scalar, enable_if_integer<Scalar> = 0>
inline Eigen::VectorX<Scalar> poly_modinv(const Eigen::VectorX<Scalar> &P,
                                          const Eigen::VectorX<Scalar> &MOD);

/* Polynomial modular inverse, could be numerically unstable */
template <typename Scalar, enable_if_not_integer<Scalar> = 0>
inline Eigen::VectorX<Scalar> poly_modinv(const Eigen::VectorX<Scalar> &P,
                                          const Eigen::VectorX<Scalar> &MOD);

/* Content of the polynomial, i.e. the gcd of all coefficients */
template <typename Scalar>
inline Scalar poly_content(const Eigen::VectorX<Scalar> &P);

/* Primitive part of a polynomial, i.e. P / poly_content(P) */
template <typename Scalar>
inline Eigen::VectorX<Scalar> poly_pp(const Eigen::VectorX<Scalar> &P);

/* Evaluate P(x) at the matrix X, i.e. P(X) */
template <typename Scalar>
inline Eigen::MatrixX<Scalar> poly_eval_mat(const Eigen::VectorX<Scalar> &P,
                                            const Eigen::MatrixX<Scalar> &X);

/* Polynomial composition P(Q(x)) (mod M(x)), to avoid exponential growth of
 * degree */
template <typename Scalar>
inline Eigen::VectorX<Scalar> poly_composition(const Eigen::VectorX<Scalar> &P,
                                               const Eigen::VectorX<Scalar> &X,
                                               const Eigen::VectorX<Scalar> &M);

template <typename Scalar>
inline Eigen::VectorX<Scalar> derivative(const Eigen::VectorX<Scalar> &P);

template <typename Scalar>
inline Eigen::VectorX<Scalar> integral(const Eigen::VectorX<Scalar> &P);

/* Compute the characteristic polynomial numerically with type 'MidType'
    => If there is only one type, it will be inferede automatically */
template <typename InType, typename MidType = InType,
          typename OutType = MidType>
inline Eigen::VectorX<OutType> charpoly(const Eigen::MatrixX<InType> &A,
                                        const bool round_flag = false);

template <typename Scalar>
inline std::string poly2string(const Eigen::VectorX<Scalar> &p);

#include "polynomial.cpp"

#endif  // POLYNOMIAL_HPP
