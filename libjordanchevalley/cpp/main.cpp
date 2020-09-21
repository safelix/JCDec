#include "io.hpp"
#include "jordanchevalley.hpp"
#include "polynomial.hpp"

#include <Eigen/Dense>
#include <iostream>  // std::cout

using namespace std;
using namespace Eigen;
namespace mp = boost::multiprecision;

/* Type Definitions */
namespace Eigen {
	
template <typename ScalarType>
using VectorX = typename Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

template <typename ScalarType>
using MatrixX =
    typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

typedef VectorX<mp::cpp_rational> VectorXq;
typedef MatrixX<mp::cpp_rational> MatrixXq;

}  // namespace Eigen

int main(int argc, char *argv[]) {
  srand((unsigned int)time(0));

  /*
  if (argc == 1)
          return 1;

  string path = string(argv[1]);
  MatrixXd A = readCSV(path);

  */

  MatrixXd A;
  /*
  // Example 1
  A = MatrixXd::Identity(4, 4);

  // Example 2
  A.diagonal(1) << 1, 1, 0;

  // Example 3
  A.resize(3, 3);
  A << 0, 0, 0, 0, 1, 1, 1, 1, 1;

  // Example 4
  A = MatrixXd::Random(3,3).array().round();

  */
  // Example 4
  A = readCSV("tests/CoutyEsterleZarouf.csv");

  JCDec jc(A, false, 2);
  cout << "timing: " << jc.get_timing() << endl;
  cout << "maxBitsize: " << jc.get_maxBitsize() << endl;

  // compute the characteristics polynomial
  jc.compute_chiA("mp::oct", true);
  cout << "timing: " << jc.get_timing() << endl;
  cout << "maxBitsize: " << jc.get_maxBitsize() << endl;

  // compute the minimal polynomial
  jc.compute_muD();
  cout << "timing: " << jc.get_timing() << endl;
  cout << "maxBitsize: " << jc.get_maxBitsize() << endl << endl;

  jc.compute_muD(true);  // with resultant
  cout << "timing: " << jc.get_timing() << endl;
  cout << "maxBitsize: " << jc.get_maxBitsize() << endl;

  // compute the inverse
  jc.compute_inv();
  cout << "timing: " << jc.get_timing() << endl;
  cout << "maxBitsize: " << jc.get_maxBitsize() << endl;

  // compute the chevalley polynomial
  jc.compute_chev();
  cout << "timing: " << jc.get_timing() << endl;
  cout << "maxBitsize: " << jc.get_maxBitsize() << endl;

  // compute the diagonalizable part
  jc.compute_D();
  cout << "timing: " << jc.get_timing() << endl;
  cout << "maxBitsize: " << jc.get_maxBitsize() << endl;

  // compute the nillpotent part
  jc.compute_N();
  cout << "timing: " << jc.get_timing() << endl;
  cout << "maxBitsize: " << jc.get_maxBitsize() << endl;

  // Print overall timings
  cout << endl << "summary timings: " << endl << jc.get_timings() << endl;
  cout << "summary maxBitsizes: " << endl
       << jc.get_maxBitsizes() << endl
       << endl;

  // check the cayleyhamilton property
  jc.check_cayleyhamilton();
  cout << "timing: " << jc.get_timing() << endl;

  // chevk the nillpotency property
  jc.check_nillpotency();
  cout << "timing: " << jc.get_timing() << endl;

  // check the commutativity property
  jc.check_commutativity();
  cout << "timing: " << jc.get_timing() << endl;

  cout << endl << ">> done <<" << endl;
  cout << endl << "A = " << endl << jc.get_A() << endl;
  cout << "chi_A(x) = " << poly2string(jc.get_chiA()) << endl;
  cout << "mu_D(x) = " << poly2string(jc.get_muD()) << endl;
  cout << "inv(x) = " << poly2string(jc.get_inv()) << endl;
  cout << "chev(x) = " << poly2string(jc.get_chev()) << endl;
  cout << "D = " << endl << jc.get_D() << endl;
  cout << "N = " << endl << jc.get_N() << endl;

  writeCSV(jc.get_N(), "tests/CoutyEsterleZarouf_res.csv");

  return 0;
}