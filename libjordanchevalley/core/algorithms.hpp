#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP

#include <Eigen/Dense>

namespace Eigen {
	
template <typename ScalarType>
using VectorX = typename Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

template <typename ScalarType>
using MatrixX =
    typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

}  // namespace Eigen

/* Remove leading zero coefficients of the polynomial p(x) */
template <typename Scalar>
inline void remove_leading_zeros(Eigen::VectorX<Scalar> &p);

/* Compute the characteristic polynomial as described in:
        La Budde's Method For Computing Characteristic Polynomials
        by Rizwana Rehman and Ilse C.F. Ipsen (2011)
        (https://arxiv.org/pdf/1104.3769.pdf) */
template <typename Scalar>
inline void labudde(const Eigen::MatrixX<Scalar> &A, Eigen::VectorX<Scalar> &p);

/* Compute the multiplication of two polynomials P(x)Q(X) = r(x)
    using the 'standard' long multiplication algorithm. */
template <typename Scalar>
inline void long_multiplication(const Eigen::VectorX<Scalar> &P,
                                const Eigen::VectorX<Scalar> &Q,
                                Eigen::VectorX<Scalar> &r);

/* Compute the euclidian division U(X) = V(X)q(X) + r(X) as described in:
        The Art of Computer Programming: Seminumerical Algorithms (Vol. 2)
        by Donald E. Knuth (1998): Chapter 4.6.1, Algorithm D */
template <typename Scalar>
inline void long_division(const Eigen::VectorX<Scalar> &U,
                          const Eigen::VectorX<Scalar> &V,
                          Eigen::VectorX<Scalar> &q, Eigen::VectorX<Scalar> &r);

/* Compute the pseudo division c*U(X) = V(X)q(X) + r(X) as described in:
        The Art of Computer Programming: Seminumerical Algorithms (Vol. 2)
        by Donald E. Knuth (1998): Chapter 4.6.1, Algorithm R */
template <typename Scalar>
inline void pseudo_division(const Eigen::VectorX<Scalar> &U,
                            const Eigen::VectorX<Scalar> &V,
                            Eigen::VectorX<Scalar> &q,
                            Eigen::VectorX<Scalar> &r);

/* Compute the greatest common divisor gcd(U(X), V(X)) as described in
        The Art of Computer Programming: Seminumerical Algorithms (Vol. 2)
        Chapter 4.5.2, Algorithm A (Modern Euclidean Algorithm)
        by Donald E. Knuth (1998) */
template <typename Scalar>
inline void euclidean_algo(const Eigen::VectorX<Scalar> &U,
                           const Eigen::VectorX<Scalar> &V,
                           Eigen::VectorX<Scalar> &gcd);

/* Compute the gcd and the bezout coefficients as descirbed in
        The Art of Computer Programming: Seminumerical Algorithms (Vol. 2)
        Chapter 4.5.2, Algorithm X (Extended Euclid's Algorithm)
        by Donald E. Knuth (1998) */
template <typename Scalar>
inline void ext_euclidean_algo(const Eigen::VectorX<Scalar> &U,
                               const Eigen::VectorX<Scalar> &V,
                               Eigen::VectorX<Scalar> &gcd,
                               Eigen::VectorX<Scalar> &u1,
                               Eigen::VectorX<Scalar> &v1);

/* Compute the primitive gcd over a unique factorization domain as described in
        The Art of Computer Programming: Seminumerical Algorithms (Vol. 2)
        Chapter 4.5.2, Algorithm E (Generalized Euclidean Algorithm)
        by Donald E. Knuth (1998) */
template <typename Scalar>
inline void euclidean_algo_ufd(const Eigen::VectorX<Scalar> &U,
                               const Eigen::VectorX<Scalar> &V,
                               Eigen::VectorX<Scalar> &gcd);

/* Compute the gcd and the bezout coefficients over a unigwue factorization
   domain. I don't know of any bounds for the bitsize coefficients of the
   results. Exponential growth could be possible! */
template <typename Scalar>
inline void ext_euclidean_algo_ufd(const Eigen::VectorX<Scalar> &U,
                                   const Eigen::VectorX<Scalar> &V,
                                   Eigen::VectorX<Scalar> &gcd,
                                   Eigen::VectorX<Scalar> &u1,
                                   Eigen::VectorX<Scalar> &v1);

/* Compute a greatest common divisor gcd(U(X), V(X)) as described in
        The Art of Computer Programming: Seminumerical Algorithms (Vol. 2)
        Chapter 4.5.2, Algorithm C (GCD over a unique factorization domain)
        by Donald E. Knuth (1998) */
template <typename Scalar>
inline void subresultant_algo(const Eigen::VectorX<Scalar> &U,
                              const Eigen::VectorX<Scalar> &V,
                              Eigen::VectorX<Scalar> &gcd);

/* Evaluate the polynomial P on the matrix X using horner method. */
template <typename Scalar>
inline void horner_method_mat(const Eigen::VectorX<Scalar> &P,
                              const Eigen::MatrixX<Scalar> &X,
                              Eigen::MatrixX<Scalar> &y);

/* Evaluate the polynomial P on the polynomial Q using the horner method mod M
 */
template <typename Scalar>
inline void horner_method_polymod(const Eigen::VectorX<Scalar> &P,
                                  const Eigen::VectorX<Scalar> &Q,
                                  const Eigen::VectorX<Scalar> &M,
                                  Eigen::VectorX<Scalar> &y);

#include "algorithms.cpp"

#endif  // ALGORITHMS_HPP
