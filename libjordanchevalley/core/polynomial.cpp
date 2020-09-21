#include "polynomial.hpp"
#include "algorithms.hpp"
#include "lowlevel.hpp"

using namespace Eigen;

/* Type Definitions */
namespace Eigen {

template <typename ScalarType>
using VectorX = typename Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

template <typename ScalarType>
using MatrixX =
    typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

template <typename ScalarType>
using ArrayX = typename Eigen::Array<ScalarType, Eigen::Dynamic, 1>;

template <typename ScalarType>
using ArrayXX =
    typename Eigen::Array<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;
	
}  // namespace Eigen

template <typename Scalar>
bool poly_isConst(const VectorX<Scalar> &P) {
  return P.size() == 1;
}

template <typename Scalar>
bool poly_isZero(const VectorX<Scalar> &P) {
  return (poly_isConst(P) && P(0) == 0);
}

template <typename Scalar>
bool poly_isOne(const VectorX<Scalar> &P) {
  return (poly_isConst(P) && P(0) == 1);
}

template <typename Scalar>
bool poly_isEqual(const VectorX<Scalar> &P, const VectorX<Scalar> &Q) {
  if (P.size() != Q.size()) return false;

  for (Index i = 0; i < P.size(); i++)
    if (P(i) != Q(i)) return false;

  return true;
}

/* Addition of polynomials, stored in a vector */
template <typename Scalar>
VectorX<Scalar> poly_add(const VectorX<Scalar> &P, const VectorX<Scalar> &Q) {
  Index n = std::max(P.size(), Q.size());

  VectorX<Scalar> res = VectorX<Scalar>::Zero(n);
  res.head(P.size()) += P;
  res.head(Q.size()) += Q;
  remove_leading_zeros(res);

  return res;
}

/* Addition of polynomials, stored in a vector */
template <typename Scalar>
VectorX<Scalar> poly_sub(const VectorX<Scalar> &P, const VectorX<Scalar> &Q) {
  Index n = std::max(P.size(), Q.size());

  VectorX<Scalar> res = VectorX<Scalar>::Zero(n);
  res.head(P.size()) += P;
  res.head(Q.size()) -= Q;
  remove_leading_zeros(res);

  return res;
}

template <typename Scalar>
VectorX<Scalar> poly_mult(const VectorX<Scalar> &P, const VectorX<Scalar> &Q) {
  VectorX<Scalar> res;
  long_multiplication(P, Q, res);

  return res;
}

template <typename Scalar, enable_if_integer<Scalar>>
VectorX<Scalar> poly_div(const VectorX<Scalar> &P, const VectorX<Scalar> &Q,
                         const bool pseudo) {
  VectorX<Scalar> res, tmp;
  if (pseudo)
    pseudo_division(P, Q, res, tmp);
  else
    pseudo_division(P, Q, res, tmp);

  return res;
}

template <typename Scalar, enable_if_not_integer<Scalar>>
VectorX<Scalar> poly_div(const VectorX<Scalar> &P, const VectorX<Scalar> &Q,
                         const bool pseudo) {
  VectorX<Scalar> res, tmp;
  if (pseudo)
    pseudo_division(P, Q, res, tmp);
  else
    long_division(P, Q, res, tmp);

  return res;
}

template <typename Scalar, enable_if_integer<Scalar>>
VectorX<Scalar> poly_mod(const VectorX<Scalar> &P, const VectorX<Scalar> &Q,
                         const bool pseudo) {
  VectorX<Scalar> res, tmp;
  if (pseudo)
    pseudo_division(P, Q, tmp, res);
  else
    pseudo_division(P, Q, tmp, res);

  return res;
}

template <typename Scalar, enable_if_not_integer<Scalar>>
VectorX<Scalar> poly_mod(const VectorX<Scalar> &P, const VectorX<Scalar> &Q,
                         const bool pseudo) {
  VectorX<Scalar> res, tmp;
  if (pseudo)
    pseudo_division(P, Q, tmp, res);
  else
    long_division(P, Q, tmp, res);

  return res;
}

template <typename Scalar, enable_if_integer<Scalar>>
VectorX<Scalar> poly_quorem(const VectorX<Scalar> &P, const VectorX<Scalar> &Q,
                            VectorX<Scalar> &rem, const bool pseudo) {
  VectorX<Scalar> quot;
  if (pseudo)
    pseudo_division(P, Q, quot, rem);
  else
    pseudo_division(P, Q, quot, rem);

  return quot;
}

template <typename Scalar, enable_if_not_integer<Scalar>>
VectorX<Scalar> poly_quorem(const VectorX<Scalar> &P, const VectorX<Scalar> &Q,
                            VectorX<Scalar> &rem, const bool pseudo) {
  VectorX<Scalar> quot;
  if (pseudo)
    pseudo_division(P, Q, quot, rem);
  else
    long_division(P, Q, quot, rem);

  return quot;
}

template <typename Scalar, enable_if_integer<Scalar>>
VectorX<Scalar> poly_gcd(const VectorX<Scalar> &P, const VectorX<Scalar> &Q) {
  VectorX<Scalar> gcd;
  subresultant_algo(P, Q, gcd);

  return gcd;
}

template <typename Scalar, enable_if_not_integer<Scalar>>
VectorX<Scalar> poly_gcd(const VectorX<Scalar> &P, const VectorX<Scalar> &Q,
                         const bool resultant) {
  VectorX<Scalar> gcd;
  if (resultant)
    subresultant_algo(P, Q, gcd);
  else
    euclidean_algo(P, Q, gcd);

  return gcd;
}

template <typename Scalar, enable_if_integer<Scalar>>
VectorX<Scalar> poly_modinv(const VectorX<Scalar> &P,
                            const VectorX<Scalar> &MOD) {
  VectorX<Scalar> gcd, tmp, inv;
  ext_euclidean_algo_ufd(P, MOD, gcd, inv, tmp);
  eigen_assert(gcd.size() == 1 && gcd(0) == 1);  // TODO: precision?

  return inv;
}

template <typename Scalar, enable_if_not_integer<Scalar>>
VectorX<Scalar> poly_modinv(const VectorX<Scalar> &P,
                            const VectorX<Scalar> &MOD) {
  VectorX<Scalar> gcd, tmp, inv;
  ext_euclidean_algo(P, MOD, gcd, inv, tmp);
  eigen_assert(gcd.size() == 1 && gcd(0) == 1);  // TODO: precision?

  return inv;
}

template <typename Scalar>
Scalar poly_content(const VectorX<Scalar> &P) {
  return P.redux([](Scalar a, Scalar b) -> Scalar { return gcd(a, b); });
}

template <typename Scalar>
VectorX<Scalar> poly_pp(const VectorX<Scalar> &P) {
  if (poly_content(P) != 0)
    return P / poly_content(P);
  else
    return P;
}

template <typename Scalar>
MatrixX<Scalar> poly_eval_mat(const VectorX<Scalar> &P,
                              const MatrixX<Scalar> &X) {
  MatrixX<Scalar> res;
  horner_method_mat(P, X, res);

  return res;
}

template <typename Scalar>
VectorX<Scalar> poly_composition(const VectorX<Scalar> &P,
                                 const VectorX<Scalar> &X,
                                 const VectorX<Scalar> &M) {
  VectorX<Scalar> res;
  horner_method_polymod(P, X, M, res);

  return res;
}

/* Derivative of a polynomial stored in a vector */
template <typename Scalar>
VectorX<Scalar> derivative(const VectorX<Scalar> &P) {
  Index n = P.size() - 1;

  ArrayX<Scalar> res = ArrayX<Scalar>::LinSpaced(n, 1, n);
  res *= P.tail(n).array();

  return res;
}

/* Integral of a polynomial stored in a vector */
template <typename Scalar>
VectorX<Scalar> integral(const VectorX<Scalar> &P) {
  Index n = P.size() + 1;

  ArrayX<Scalar> res(n);
  res << 0, P.array();
  res.tail(n - 1) /= ArrayX<Scalar>::LinSpaced(n - 1, 1, n - 1);
  return res;
}

/* Compute the characteristic polynomial */
template <typename InType, typename MidType, typename OutType>
VectorX<OutType> charpoly(const MatrixX<InType> &A, const bool round_flag) {
  eigen_assert(A.rows() == A.cols());
  Index n = A.rows();

  VectorX<MidType> p(n);
  labudde(A.template cast<MidType>().eval(), p);

  if (round_flag) p = p.unaryExpr([](MidType p_i){return round(p_i); });

  return p.template cast<OutType>();
}

template <typename Scalar>
std::string poly2string(const VectorX<Scalar> &p) {
  std::ostringstream stream;
  Index deg = p.size() - 1;
  for (Index i = deg; i > -1; i--) {
    stream << p(i);
    if (i == 0)
      stream << "";
    else if (i == 1)
      stream << "x + ";
    else {
      stream << "x";
      for (char &c : std::to_string(i)) {
        if (c == '0') stream << "⁰";
        if (c == '1') stream << "¹";
        if (c == '2') stream << "²";
        if (c == '3') stream << "³";
        if (c == '4') stream << "⁴";
        if (c == '5') stream << "⁵";
        if (c == '6') stream << "⁶";
        if (c == '7') stream << "⁷";
        if (c == '8') stream << "⁸";
        if (c == '9') stream << "⁹";
      }
      stream << " + ";
    }
  }

  return stream.str();
}
