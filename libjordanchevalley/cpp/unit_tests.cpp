#include "algorithms.hpp"
#include "polynomial.hpp"

#include <Eigen/Dense>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/eigen.hpp>

#include <iostream>  // std::cout

using namespace std;
using namespace Eigen;
namespace mp = boost::multiprecision;

namespace Eigen {

template <typename ScalarType>
using VectorX = typename Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

template <typename ScalarType>
using MatrixX =
    typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

typedef VectorX<mp::cpp_rational> VectorXq;
typedef MatrixX<mp::cpp_rational> MatrixXq;

}  // namespace Eigen

template <typename Scalar>
bool test_addsub(const VectorX<Scalar> &P, const VectorX<Scalar> &Q) {
  return poly_isEqual(P, poly_sub(poly_add(P, Q), Q));
}

template <typename Scalar>
bool test_multdiv(const VectorX<Scalar> &P, const VectorX<Scalar> &Q) {
  return poly_isEqual(P, poly_div(poly_mult(P, Q), Q));
}

template <typename Scalar>
bool test_divmult(const VectorX<Scalar> &P, const VectorX<Scalar> &Q,
                  bool pseudo) {
  VectorX<Scalar> quot, rem, check;

  quot = poly_quorem(P, Q, rem);
  check = poly_add(poly_mult(quot, Q), rem);

  if (pseudo || std::numeric_limits<Scalar>::is_integer)
    return poly_isConst(poly_div(check, P)) && poly_isZero(poly_mod(check, P));
  else
    return poly_isEqual(P, check);
}

template <typename Scalar>
bool test_derivint(const VectorX<Scalar> &P) {
  return poly_isConst(poly_sub(P, integral(derivative(P))));
}

template <typename Scalar>
bool test_gcd(const VectorX<Scalar> &P, const VectorX<Scalar> &Q,
              bool resultant = false) {
  if (!poly_isZero(poly_mod(P, poly_gcd(P, Q, resultant)))) return false;

  if (!poly_isZero(poly_mod(Q, poly_gcd(P, Q, resultant)))) return false;

  return true;
}

template <typename Scalar>
bool test_modinv(const VectorX<Scalar> &P, const VectorX<Scalar> &Q) {
  VectorX<Scalar> invP = poly_modinv(P, Q);
  return poly_isOne(poly_mod(poly_mult(invP, P), Q));
}

template <typename Scalar>
bool test_charpoly(const MatrixX<Scalar> &A) {
  VectorX<Scalar> chi_A = charpoly(A);
  return poly_eval_mat(chi_A, A).isZero();
}

int main() {
  // https://en.wikipedia.org/wiki/Polynomial_long_division#Polynomial_long_division
  cout << endl;
  VectorXq p, q;
  p.resize(4);
  q.resize(2);
  p << -4, 0, -2, 1;
  q << -9, 3;
  cout << "p(x) = " << poly2string(p) << endl;
  cout << "q(x) = " << poly2string(q) << endl;

  cout << "long multiplication: prod(x) = p(x) * q(x) " << flush;
  if (test_multdiv(p, q))
    cout << ">> passed" << endl;
  else {
    cout << ">> failed" << endl;
    cout << ">> prod(x) = " << poly2string(poly_mult(p, q)) << endl;
  }

  cout << "long division: p(x) = q(x)*quot(x) + rem(x) " << flush;
  if (test_divmult(p, q, false))
    cout << ">> passed" << endl;
  else {
    cout << ">> failed " << endl;
    cout << ">> quot(x) = " << poly2string(poly_div(p, q)) << endl;
    cout << ">> rem(x) = " << poly2string(poly_mod(p, q)) << endl;
  }

  cout << "pseudo division: const*p(x) = q(x)*quot(x) + rem(x) " << flush;
  if (test_divmult(p, q, true))
    cout << ">> passed" << endl;
  else {
    cout << ">> failed" << endl;
    cout << ">> quot(x) = " << poly2string(poly_div(p, q, true)) << endl;
    cout << ">> rem(x) = " << poly2string(poly_mod(p, q, true)) << endl;
  }

  cout << "integral of derivative: int_0^x( d/dx( p(x) ) ) " << flush;
  if (test_derivint(p))
    cout << ">> passed" << endl;
  else {
    cout << ">> failed" << endl;
    cout << ">> p'(x) ) = " << poly2string(derivative(p)) << endl;
    cout << ">> int_0^x( p'(x) ) =  " << poly2string(integral(derivative(p)))
         << endl;
  }

  // https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Pseudo-remainder_sequences
  cout << endl;
  /*
  p.resize(7);
  q.resize(9);
  p << 21, -9, -4, 0, 5, 0, 3;
  q << -5, 2, 8, -3, -3, 0, 1, 0, 1;
  */
  p.resize(1);
  q.resize(2);
  p << 1;
  q << 1, 1;
  cout << "p(x) = " << poly2string(p) << endl;
  cout << "q(x) = " << poly2string(q) << endl;

  cout << "greatest common divisor gcd(x) = gcd(p(x), q(x)) " << flush;
  if (test_gcd(p, q))
    cout << ">> passed" << endl;
  else {
    cout << ">> failed" << endl;
    cout << ">> gcd(x) = " << poly2string(poly_gcd(p, q)) << endl;
    cout << ">> p(x) % gcd(x) = " << poly2string(poly_mod(p, poly_gcd(p, q)))
         << endl;
    cout << ">> q(x) % gcd(x) = " << poly2string(poly_mod(q, poly_gcd(p, q)))
         << endl;
  }

  cout << "greatest common divisor gcd(x) = gcd(p(x), q(x))   (resultant algo)"
       << flush;
  if (test_gcd(p, q, true))
    cout << ">> passed" << endl;
  else {
    cout << ">> failed" << endl;
    cout << ">> gcd(x) = " << poly2string(poly_gcd(p, q, true)) << endl;
    cout << ">> p(x) % gcd(x) = "
         << poly2string(poly_mod(p, poly_gcd(p, q, true))) << endl;
    cout << ">> q(x) % gcd(x) = "
         << poly2string(poly_mod(q, poly_gcd(p, q, true))) << endl;
  }

  cout << "modular inverse inv(x) * p(x) % q(x) = 1 " << flush;
  if (test_modinv(p, q))
    cout << ">> passed" << endl;
  else {
    cout << ">> failed" << endl;
    cout << ">> inv(x) = " << poly2string(poly_modinv(p, q)) << endl;
    cout << ">> inv(x) * p(x) % q(x) = "
         << poly2string(poly_mod(poly_mult(poly_modinv(p, q), p), q)) << endl;
  }

  /* https://en.wikipedia.org/wiki/Characteristic_polynomial */
  cout << endl;
  MatrixX<mp::cpp_bin_float_100> A(2, 2);
  A << 2, 1, -1, 0;
  cout << "A = " << endl << A << endl;

  cout << "characteristic polynomial: chi_A(A) = 0 " << flush;
  if (test_charpoly(A))
    cout << ">> passed" << endl;
  else {
    cout << ">> failed" << endl;
    cout << "chi_A(x) = " << poly2string(charpoly(A)) << endl;
    cout << "chi_A(A) = " << endl << poly_eval_mat(charpoly(A), A) << endl;
  }

  return 0;
}