#include "jordanchevalley.hpp"

#include <Eigen/Dense>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <chrono>
#include <iostream>  // std::cout

#include "polynomial.hpp"

using namespace Eigen;
namespace mp = boost::multiprecision;

namespace boost {
namespace multiprecision {

typedef number<backends::cpp_bin_float<1000, backends::digit_base_10,
                                       std::allocator<void>>>
    cpp_bin_float_1000;
}

}  // namespace boost

/*
--------------------------------- Constructors ---------------------------------
*/

JCDec::JCDec(const MatrixXd& mat, bool rowmajor, int verbosity) {
  // start time
  auto start_time = std::chrono::steady_clock::now();

  // output
  if (1 <= verbosity)
    std::cout << "Compute the Jordan-Chevalley decomposition of A" << std::endl;

  // store matrix into A
  if (rowmajor)
    this->A = mat.transpose().cast<mp::cpp_rational>().eval();
  else
    this->A = mat.cast<mp::cpp_rational>().eval();

  // initiallize other attributes
  this->size = static_cast<int>(mat.rows());
  this->no_iter = -1;
  this->verbosity = verbosity;

  // output
  if (3 <= verbosity)
    std::cout << ">> A = " << std::endl << this->A << std::endl;

  // enforce squareness
  if (mat.rows() != mat.cols())
    throw std::invalid_argument("JCDec(A): A MUST BE SQUARE");

  // stop timing
  auto stop_time = std::chrono::steady_clock::now();
  this->timing = std::chrono::duration<double>(stop_time - start_time).count();
  JCDec::store_timing(0, this->timing);
  last_step = "A";
}

JCDec::JCDec(const Map<MatrixXd>& map, bool rowmajor, int verbosity) {
  // start time
  auto start_time = std::chrono::steady_clock::now();

  // output
  if (1 <= verbosity)
    std::cout << "Compute the Jordan-Chevalley decomposition of A" << std::endl;

  // store mapped array into A
  if (rowmajor)
    this->A = map.transpose().cast<mp::cpp_rational>();
  else
    this->A = map.cast<mp::cpp_rational>();

  // initiallize other attributes
  this->size = static_cast<int>(map.rows());
  this->no_iter = -1;
  this->verbosity = verbosity;

  // output
  if (3 <= verbosity)
    std::cout << ">> A = " << std::endl << this->A << std::endl;

  // enforce squarness
  if (map.rows() != map.cols())
    throw std::invalid_argument("JCDec(A): A MUST BE SQUARE");

  // stop timing
  auto stop_time = std::chrono::steady_clock::now();
  this->timing = std::chrono::duration<double>(stop_time - start_time).count();
  JCDec::store_timing(0, this->timing);
  last_step = "A";
}

/*
-------------------------------- Compute methods -------------------------------
*/

template <typename precision>
JCDec& JCDec::compute_chiA(bool round) {
  // start timing
  auto start_time = std::chrono::steady_clock::now();

  // output
  if (1 <= verbosity)
    std::cout << std::endl
              << "1. compute the characteristic polynomial of A, chi_A(x) = "
                 "PROD (x - lambda_i)^(a_i)"
              << std::endl;

  // actual computation
  this->chiA =
      charpoly<mp::cpp_rational, precision, mp::cpp_rational>(this->A, round);

  // output
  if (2 <= verbosity)
    std::cout << ">> chi_A(x) = " << poly2string(this->chiA) << std::endl;

  // stop timing
  auto stop_time = std::chrono::steady_clock::now();
  this->timing = std::chrono::duration<double>(stop_time - start_time).count();
  JCDec::store_timing(1, this->timing);
  last_step = "chiA";

  return *this;
}

JCDec& JCDec::compute_chiA(std::string precision, bool round) {
  /* hardaware floating-point types */
  if (precision == "float") return JCDec::compute_chiA<float>(round);

  if (precision == "double") return JCDec::compute_chiA<double>(round);

  if (precision == "long double")
    return JCDec::compute_chiA<long double>(round);

  /* software floating-point types */
  if (precision == "mp::single")
    return JCDec::compute_chiA<mp::cpp_bin_float_single>(round);

  if (precision == "mp::double")
    return JCDec::compute_chiA<mp::cpp_bin_float_double>(round);

  if (precision == "mp::quad")
    return JCDec::compute_chiA<mp::cpp_bin_float_quad>(round);

  if (precision == "mp::oct")
    return JCDec::compute_chiA<mp::cpp_bin_float_oct>(round);

  if (precision == "mp::50")
    return JCDec::compute_chiA<mp::cpp_bin_float_50>(round);

  if (precision == "mp::100")
    return JCDec::compute_chiA<mp::cpp_bin_float_100>(round);

  if (precision == "mp::1000")
    return JCDec::compute_chiA<mp::cpp_bin_float_1000>(round);

  throw std::invalid_argument("compute_chi_A(): unkown precision specifier");
  return *this;
}

JCDec& JCDec::compute_muD(bool resultant) {
  if (this->chiA.size() == 0) JCDec::compute_chiA("mp::oct", true);

  // start timing
  auto start_time = std::chrono::steady_clock::now();

  // Output
  if (1 <= verbosity)
    std::cout << std::endl
              << "2. compute mu_D(x) = chi_A(x) / gcd( chi_A(x), chi_A'(x) ) "
                 "=> mu_D(x) = PROD (x - lambda_i)"
              << std::endl;

  // actual computation
  VectorXq rem;
  this->muD = poly_quorem(
      this->chiA, poly_gcd(this->chiA, derivative(this->chiA), resultant), rem);

  // output
  if (2 <= this->verbosity)
    std::cout << ">> mu_D(x) = " << poly2string(this->muD) << std::endl;
  if (3 <= this->verbosity)
    std::cout << ">> remander of division = " << poly2string(rem) << std::endl;

  // stop timing
  auto stop_time = std::chrono::steady_clock::now();
  this->timing = std::chrono::duration<double>(stop_time - start_time).count();
  JCDec::store_timing(2, this->timing);
  last_step = "muD";

  return *this;
}

JCDec& JCDec::compute_inv() {
  if (this->muD.size() == 0) JCDec::compute_muD();

  // start timing
  auto start_time = std::chrono::steady_clock::now();

  // output
  if (1 <= verbosity)
    std::cout << std::endl
              << "3. compute inv(x), s.t. muD'(x) * inv(x) = 1  (mod muD(x))"
              << std::endl;

  // actual computation
  VectorXq muDp = derivative(this->muD);
  this->inv = poly_modinv(muDp, this->muD);

  // output
  if (2 <= this->verbosity)
    std::cout << ">> inv(x) = " << poly2string(this->inv) << std::endl;
  if (3 <= this->verbosity)
    std::cout << ">> muD'(x) * inv(x) % muD(x) = "
              << poly2string(poly_mod(poly_mult(muDp, inv), muD)) << std::endl;

  // stop timing
  auto stop_time = std::chrono::steady_clock::now();
  this->timing = std::chrono::duration<double>(stop_time - start_time).count();
  JCDec::store_timing(3, this->timing);
  last_step = "inv";

  return *this;
}

JCDec& JCDec::compute_chev() {
  if (this->inv.size() == 0) JCDec::compute_inv();

  // start timing
  auto start_time = std::chrono::steady_clock::now();

  // output
  if (1 <= verbosity)
    std::cout << std::endl
              << "4. iterate untill convergence: chev(x) <- chev(x) - "
                 "muD(chev(x)) * inv(chev(x))  (mod chi_A(x))"
              << std::endl;

  // actual computation
  this->chev = (VectorXq(2) << 0, 1).finished();  // chev(x) = x
  VectorXq frac =
      poly_mod(poly_mult(this->muD, this->inv),
               this->chiA);  // frac(x) = muD(x) * inv(x) % mod chi_A(x)

  // output
  if (3 <= verbosity) {
    std::cout << ">> chev(x) = " << poly2string(this->chev) << std::endl;
    std::cout << ">> muD(x) * inv(x) = " << poly2string(frac) << std::endl;
  }

  // the chevalley iteration
  VectorXq step;
  Index l, n = this->size;
  for (l = 0; l < n; l++) {
    if (3 <= verbosity)
      std::cout << std::endl << ">> Iteration: " << l << std::endl;

    // actual computation: step(x) = frac(chev(x)) (mod chi_A(x))
    step = poly_composition(frac, this->chev, this->chiA);

    // output
    if (3 <= verbosity)
      std::cout << ">> step(x) = " << poly2string(step) << std::endl;
    if (4 <= verbosity)
      std::cout << ">> step(A) = " << std::endl
                << poly_eval_mat(step, this->A) << std::endl;

    // actual computation: chev(x) <- chev(x) - step(x)
    if (poly_isZero(step))
      break;
    else
      this->chev = poly_sub(this->chev, step);

    // output
    if (3 <= verbosity)
      std::cout << ">> chev(x) = " << poly2string(this->chev) << std::endl;
    if (4 <= verbosity)
      std::cout << ">> chev(A) = " << std::endl
                << poly_eval_mat(this->chev, this->A) << std::endl;
  }

  // store number of iteration in JCDec object
  this->no_iter = l;

  // output
  if (1 <= verbosity) std::cout << ">> no_iter: " << this->no_iter << std::endl;
  if (2 <= verbosity)
    std::cout << ">> chev(x) = " << poly2string(this->chev) << std::endl;

  // stop timing
  auto stop_time = std::chrono::steady_clock::now();
  this->timing = std::chrono::duration<double>(stop_time - start_time).count();
  JCDec::store_timing(4, this->timing);
  last_step = "chev";

  return *this;
}

JCDec& JCDec::compute_D_poly() {
  if (JCDec::chev.size() == 0) JCDec::compute_chev();

  // start timing
  auto start_time = std::chrono::steady_clock::now();

  // output
  if (1 <= verbosity)
    std::cout << std::endl
              << "5. compute diagonalizable part of A: D = chev(A)"
              << std::endl;

  // actual computation
  this->D = poly_eval_mat(this->chev, this->A);

  // output
  if (3 <= verbosity)
    std::cout << ">> D = " << std::endl << this->D << std::endl;

  // stop timing
  auto stop_time = std::chrono::steady_clock::now();
  this->timing = std::chrono::duration<double>(stop_time - start_time).count();
  JCDec::store_timing(5, this->timing);
  last_step = "D";

  return *this;
}

JCDec& JCDec::compute_D_mat() {
  if (this->muD.size() == 0) JCDec::compute_muD();

  // start timing
  auto start_time = std::chrono::steady_clock::now();

  // output
  if (1 <= verbosity)
    std::cout << std::endl
              << "5. iterate untill convergence: S <- S - muD(S) * muD'(S)^{-1}"
              << std::endl;

  // actual computation
  MatrixXq S = this->A;  // start with A
  VectorXq muDp = derivative(this->muD);

  MatrixXq step;
  Index l, n = this->size;
  for (l = 0; l < n; l++) {
    if (3 <= verbosity)
      std::cout << std::endl << ">> Iteration: " << l << std::endl;

    // actual computation: step = muD(S) * muD'(S)^{-1}
    step = poly_eval_mat(muD, S);

    // output
    if (4 <= verbosity)
      std::cout << ">> muD(S) * muD'(S)^{-1} = " << std::endl
                << step << std::endl;

    // actual computation: S <- S - step
    if (step.isZero())
      break;
    else
      S = S - step * poly_eval_mat(muDp, S).inverse();

    // output
    if (4 <= verbosity) std::cout << ">> S = " << std::endl << S << std::endl;
  }

  // store number of iteration in JCDec object
  this->no_iter = l;
  this->D = S;

  if (1 <= verbosity) std::cout << ">> no_iter: " << this->no_iter << std::endl;
  if (3 <= verbosity)
    std::cout << ">> D = " << std::endl << this->D << std::endl;

  // stop timing
  auto stop_time = std::chrono::steady_clock::now();
  this->timing = std::chrono::duration<double>(stop_time - start_time).count();
  JCDec::store_timing(3, std::numeric_limits<double>::quiet_NaN());
  JCDec::store_timing(4, std::numeric_limits<double>::quiet_NaN());
  JCDec::store_timing(5, this->timing);
  last_step = "D";


  return *this;
}

JCDec& JCDec::compute_D(bool mat) {
  if (mat)
    return compute_D_mat();
  else
    return compute_D_poly();
}

JCDec& JCDec::compute_N() {
  if (JCDec::D.size() == 0) JCDec::compute_D();

  // start timing
  auto start_time = std::chrono::steady_clock::now();

  // output
  if (1 <= verbosity)
    std::cout << std::endl
              << "6. compute nillpotent part of A: N = A - D" << std::endl;

  // actual computation
  this->N = this->A - this->D;

  // output
  if (3 <= verbosity)
    std::cout << ">> N = " << std::endl << this->N << std::endl;

  // stop timing
  auto stop_time = std::chrono::steady_clock::now();
  this->timing = std::chrono::duration<double>(stop_time - start_time).count();
  JCDec::store_timing(6, this->timing);
  last_step = "N";

  return *this;
}

JCDec& JCDec::compute() { return compute_N(); }

/*
--------------------------------- Check methods --------------------------------
*/

bool JCDec::is_trivial() {
  timing = 0;

  if (this->N.size() == 0)
    throw std::invalid_argument("is_trivial(): please call compute() first");

  return this->N.isZero();
}

double JCDec::check_cayleyhamilton() {
  if (this->chiA.size() == 0)
    throw std::invalid_argument(
        "check_cayleyhamilton(): chi_A(x) does not exist");

  // start timing
  auto start_time = std::chrono::steady_clock::now();

  // output
  if (1 <= verbosity)
    std::cout
        << "Check whether the Cayley-Hamilton property holds: chi_A(A) = 0"
        << std::endl;

  // actual computation
  mp::cpp_rational err = poly_eval_mat(this->chiA, A).squaredNorm();

  // output
  if (2 <= verbosity)
    std::cout << ">> ( || chi_A(A) ||_F )^2 = " << err << std::endl;

  // stop timing
  auto stop_time = std::chrono::steady_clock::now();
  this->timing = std::chrono::duration<double>(stop_time - start_time).count();

  return static_cast<double>(err);
}

double JCDec::check_nillpotency() {
  if (this->N.size() == 0)
    throw std::invalid_argument("check_nipplotency(): N does not exist");

  // start timing
  auto start_time = std::chrono::steady_clock::now();

  // output
  if (1 <= verbosity)
    std::cout << "Check whether N is nillpotent: N^n = 0"
              << std::endl;

  // actual computation
  mp::cpp_rational err = pow(N, this->size).squaredNorm();

  // output
  if (2 <= verbosity)
    std::cout << ">> ( || N^n ||_F )^2 = " << err << std::endl;

  // stop timing
  auto stop_time = std::chrono::steady_clock::now();
  this->timing = std::chrono::duration<double>(stop_time - start_time).count();

  return static_cast<double>(err);
}

double JCDec::check_commutativity() {
  if (this->N.size() == 0 || this->D.size() == 0)
    throw std::invalid_argument("check_commutativity(): D or N does not exist");

  // start timing
  auto start_time = std::chrono::steady_clock::now();

  // output
  if (1 <= verbosity)
    std::cout << "Check whether D and N commute: D*N - N*D = 0" << std::endl;

  // actual computation
  mp::cpp_rational err = (D * N - N * D).squaredNorm();

  // output
  if (2 <= verbosity)
    std::cout << ">> ( D*N - N*D ||_F )^2 = " << err << std::endl;

  // stop timing
  auto stop_time = std::chrono::steady_clock::now();
  this->timing = std::chrono::duration<double>(stop_time - start_time).count();

  return static_cast<double>(err);
}

/*
-------------------------------- Helper functions
-------------------------------
*/

void JCDec::store_timing(int pos, double timing) {
  this->timings.conservativeResize(pos + 1);
  this->timings(pos) = timing;
}

double JCDec::get_timing(std::string step_name) {
  if (step_name == "") return this->timing;

  if (step_name == "A") return this->timings(0);

  if (step_name == "chiA") return this->timings(1);

  if (step_name == "muD") return this->timings(2);

  if (step_name == "inv") return this->timings(3);

  if (step_name == "chev") return this->timings(4);

  if (step_name == "D") return this->timings(5);

  if (step_name == "N") return this->timings(6);

  throw std::invalid_argument("get_timing(): unkown step_name");
}

RowVectorXd JCDec::get_timings(int steps) {
  if (steps < 1)
    throw std::invalid_argument(
        "get_timings(): please specify more than 0 steps");

  steps = std::min(steps, static_cast<int>(timings.size()));
  return timings.head(steps);
}

inline double bitsize(mp::cpp_rational frac) {
  double num = static_cast<double>(mp::numerator(frac));
  double denom = static_cast<double>(mp::denominator(frac));

  return ceil(log2(abs(num) + 1)) + ceil(log2(abs(denom)));
}

double JCDec::get_maxBitsize(std::string step_name) {
  MatrixXq attr;
  if (step_name == "")
    return get_maxBitsize(last_step);

  else if (step_name == "A")
    attr = this->A;

  else if (step_name == "chiA")
    attr = this->chiA;

  else if (step_name == "muD")
    attr = this->muD;

  else if (step_name == "inv")
    attr = this->inv;

  else if (step_name == "chev")
    attr = this->chev;

  else if (step_name == "D")
    attr = this->D;

  else if (step_name == "N")
    attr = this->N;

  else
    throw std::invalid_argument("get_maxBitsize(): unkown name");

  if (attr.size() == 0)
    return std::numeric_limits<double>::quiet_NaN();
    /*
    throw std::invalid_argument("get_maxBitsize(): " + step_name +
                                " has not been computed");
    */ 

  return attr.unaryExpr(std::ref(bitsize)).maxCoeff();
}

RowVectorXd JCDec::get_maxBitsizes(int steps) {
  if (steps < 1)
    throw std::invalid_argument(
        "get_maxBitsizes(): please specify more than 0 steps");

  steps = std::min(steps, static_cast<int>(timings.size()));

  RowVectorXd maxBitsizes(steps);
  for (int i = 0; i < steps; i++)
    maxBitsizes(i) = get_maxBitsize(this->names[i]);

  return maxBitsizes;
}
