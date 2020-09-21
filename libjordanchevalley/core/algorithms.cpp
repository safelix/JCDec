#include "algorithms.hpp"

#include "lowlevel.hpp"
#include "polynomial.hpp"

using namespace Eigen;

template <typename Scalar>
void remove_leading_zeros(VectorX<Scalar> &p) {
  Index i = p.size() - 1;
  // TODO: chose adequate tolerance
  while (i > 0 && p(i) == 0.0) i--;
  p = p.head(i + 1).eval();
}

template <typename Scalar>
void labudde(const MatrixX<Scalar> &A, VectorX<Scalar> &p) {
  eigen_assert(A.rows() == A.cols());
  Index n = A.rows();

  HessenbergDecomposition<MatrixX<Scalar>> hess(A);
  MatrixX<Scalar> H = hess.matrixH();
  VectorX<Scalar> alpha = H.diagonal(0);
  VectorX<Scalar> betas = H.diagonal(-1);

  /* Let p_i(X) be the characteristic polynomial of the
    leading principal submatrix H_i of order i. The ith
    column of c stores the coefficients of p_i(X). */
  MatrixX<Scalar> c = MatrixX<Scalar>::Zero(n + 1, n + 1);

  // p_0(X) = x - alpha_0
  c(0, 0) = 1;
  for (Index i = 1; i < n + 1; i++) {
    // p_i(X) = p_{i-1}(X) * X
    c.col(i).segment(1, i) = c.col(i - 1).segment(0, i);

    // p_i(X) -= p_{i-1}(X) * alpha_i
    c.col(i) -= c.col(i - 1) * alpha(i - 1);

    // p_i(X) -= linearCombination(p_0(X), ..., p_{i-2}(X))
    if (1 < i) betas.head(i - 2) *= betas(i - 2);
    c.col(i) -= c.leftCols(i - 1) *
                betas.head(i - 1).cwiseProduct(H.col(i - 1).head(i - 1));
  }

  p = c.col(n);
}

template <typename Scalar>
void long_multiplication(const VectorX<Scalar> &P, const VectorX<Scalar> &Q,
                         VectorX<Scalar> &r) {
  Index nP = P.size();
  Index nQ = Q.size();
  Index n = nP + nQ - 1;

  r = VectorX<Scalar>::Zero(n);
  for (Index i = 0; i < nP; i++) {
    r.segment(i, nQ) += P(i) * Q;
  }
}

template <typename Scalar>
void long_division(const VectorX<Scalar> &U, const VectorX<Scalar> &V,
                   VectorX<Scalar> &q, VectorX<Scalar> &r) {
  Index degU = U.size() - 1;
  Index degV = V.size() - 1;

  r = U;
  if (degU < degV) {
    q = VectorX<Scalar>::Zero(1);
    return;
  }

  Index degq = degU - degV;
  q = VectorX<Scalar>::Zero(degq + 1);
  for (Index k = 0; k < degq + 1; k++) {
    // coefficient
    q(degq - k) = r(degU - k) / V(degV);

    // r(degq - k) => 0
    r.segment(degq - k, degV + 1) -= q(degq - k) * V;
  }

  remove_leading_zeros(r);
}

template <typename Scalar>
void pseudo_division(const VectorX<Scalar> &U, const VectorX<Scalar> &V,
                     VectorX<Scalar> &q, VectorX<Scalar> &r) {
  Index degU = U.size() - 1;
  Index degV = V.size() - 1;

  r = U;
  if (degU < degV) {
    q = VectorX<Scalar>::Zero(1);
    return;
  }

  // r *= pow(V(degV), degU - degV + 1);

  Index degq = degU - degV;
  q = VectorX<Scalar>::Zero(degq + 1);
  for (Index k = 0; k < degq + 1; k++) {
    // r(x) : v(x) = q(x) |*lc(v(x))=V(degV)
    q *= V(degV);
    r *= V(degV);

    // V(degV) divides all elements of r => exact
    q(degq - k) = r(degU - k) / V(degV);

    // r(degq - k) => 0
    r.segment(degq - k, degV + 1) -= q(degq - k) * V;
  }

  /* Optimization: do not call poly_content()
- poly_content is slow and coefficient growth moderate => faster to leave out
- Subresultant needs the pseudo-coefficient to be pow(V(degV), degU - degV +
1) Scalar cont = gcd(poly_content(q), poly_content(r)); r /= cont; q /= cont;
*/

  remove_leading_zeros(r);
}

template <typename Scalar>
void euclidean_algo(const VectorX<Scalar> &U, const VectorX<Scalar> &V,
                    VectorX<Scalar> &gcd) {
  VectorX<Scalar> u, v, r;

  u = U;
  v = V;
  while (!poly_isZero(v)) {
    r = poly_mod(u, v);
    u = v;
    v = r;
  }
  gcd = u / u(u.size() - 1);
}

template <typename Scalar>
void ext_euclidean_algo(const VectorX<Scalar> &U, const VectorX<Scalar> &V,
                        VectorX<Scalar> &gcd, VectorX<Scalar> &u1,
                        VectorX<Scalar> &v1) {
  VectorX<Scalar> u2, v2, v3, u3, q, r;
  u1 = VectorX<Scalar>::Ones(1);
  u2 = VectorX<Scalar>::Zero(1);
  u3 = U;

  v1 = VectorX<Scalar>::Zero(1);
  v2 = VectorX<Scalar>::Ones(1);
  v3 = V;

  while (!poly_isZero(v3)) {
    q = poly_quorem(u3, v3, r);
    u3 = v3;
    v3 = r;

    // temp = u2; u2 = u1 - q*u2; u1 = temp;
    u1.swap(u2);
    u2 = poly_sub(u2, poly_mult(u1, q));

    // temp = v2; v2 = v1 - q*v2; u1 = temp;
    v1.swap(v2);
    v2 = poly_sub(v2, poly_mult(v1, q));
  }

  // normalize to get monomial gcd
  u1 /= u3(u3.size() - 1);  // TODO: remove leading zeros?
  v1 /= u3(u3.size() - 1);  // TODO: remove leading zeros?
  u3 /= u3(u3.size() - 1);

  gcd = u3;
}

template <typename Scalar>
void euclidean_algo_ufd(const VectorX<Scalar> &U, const VectorX<Scalar> &V,
                        VectorX<Scalar> &gcd) {
  VectorX<Scalar> u, v, r;

  u = poly_pp(U);
  v = poly_pp(V);
  while (!poly_isZero(v)) {
    r = poly_mod(u, v);
    u = v;
    v = poly_pp(r);
  }
  gcd = poly_pp(u);
}

template <typename Scalar>
void ext_euclidean_algo_ufd(const VectorX<Scalar> &U, const VectorX<Scalar> &V,
                            VectorX<Scalar> &u3, VectorX<Scalar> &u1,
                            VectorX<Scalar> &v1) {
  Scalar cont;
  VectorX<Scalar> u2, v2, v3, q, r;
  u1 = VectorX<Scalar>::Ones(1);
  u2 = VectorX<Scalar>::Zero(1);
  u3 = U;

  v1 = VectorX<Scalar>::Zero(1);
  v2 = VectorX<Scalar>::Ones(1);
  v3 = V;

  while (!poly_isZero(v3)) {
    q = poly_quorem(u3, v3, r);
    u3 = v3;
    v3 = r;

    // temp = u2; u2 = u1 - q*u2; u1 = temp;
    u1.swap(u2);
    u2 = poly_sub(u2, poly_mult(u1, q));

    // temp = v2; v2 = v1 - q*v2; u1 = temp;
    v1.swap(v2);
    v2 = poly_sub(v2, poly_mult(v1, q));

    cont = poly_content(v3);
    cont = gcd(cont, poly_content(v2));
    cont = gcd(cont, poly_content(u2));
    if (cont != 0) {
      v3 /= cont;
      u2 /= cont;
      v2 /= cont;
    }
  }

  cont = poly_content(u3);
  cont = gcd(cont, poly_content(u1));
  cont = gcd(cont, poly_content(v1));

  u3 /= cont;
  u1 /= cont;
  v1 /= cont;
}

template <typename Scalar>
void subresultant_algo(const VectorX<Scalar> &U, const VectorX<Scalar> &V,
                       VectorX<Scalar> &gcd) {
  VectorX<Scalar> u, v, r;

  // Reduce to primitive
  u = poly_pp(U);
  v = poly_pp(V);

  // Invariant: u.size() >= v.size()
  if (u.size() < v.size()) u.swap(v);

  Scalar g = 1, h = 1;
  while (!poly_isZero(v)) {
    Index delta = u.size() - v.size();  // Invariant: delta >= 0
    eigen_assert(delta >= 0 && "subresultant: violated invariant");

    r = poly_mod(u, v, true);
    u = v;
    v = r / (g * pow(h, delta));

    g = u(u.size() - 1);
    h = (h * pow(g, delta)) / pow(h, delta);  // h = h^{ 1 - delta} * g^{delta}
  }

  gcd = poly_pp(u);
}

template <typename Scalar>
void horner_method_mat(const VectorX<Scalar> &P, const MatrixX<Scalar> &X,
                       MatrixX<Scalar> &y) {
  eigen_assert(X.rows() == X.cols());

  Index n = X.cols();
  Index deg = P.size() - 1;
  y = MatrixX<Scalar>::Zero(n, n);

  MatrixX<Scalar> id_mat = MatrixX<Scalar>::Identity(n, n);
  for (Index i = deg; i > -1; i--) y = y * X + P(i) * id_mat;
}

template <typename Scalar>
void horner_method_polymod(const VectorX<Scalar> &P, const VectorX<Scalar> &X,
                           const VectorX<Scalar> &M, VectorX<Scalar> &y) {
  Index degP = P.size() - 1;
  y = VectorX<Scalar>::Zero(1);

  for (Index i = degP; i > -1; i--) {
    y = poly_mult(y, X);
    y(0) += P(i);  // TODO: remove eval
    y = poly_mod(y, M);
  }

  remove_leading_zeros(y);
}
