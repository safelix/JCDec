#ifndef JORDANCHEVALLEY_HPP
#define JORDANCHEVALLEY_HPP

#include <Eigen/Dense>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/eigen.hpp>

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


/* Jordan-Chevalley Decomposition Module */
class JCDec {
 private:
  int size, no_iter;
  Eigen::MatrixX<mp::cpp_rational> A, D, N;
  Eigen::VectorX<mp::cpp_rational> chiA, muD, inv, chev;

  // Meta attributes
  std::string last_step;
  std::string names[7] = {"A", "chiA", "muD", "inv", "chev", "D", "N"};

  double timing;
  Eigen::RowVectorXd timings;
  inline void store_timing(int pos, double timing);

 public:
  int verbosity;

  // Constructors
  JCDec(const Eigen::MatrixXd& A, bool rowmajor = false, int verbosity = 0);
  JCDec(const Eigen::Map<Eigen::MatrixXd>& A, bool rowmajor = false,
        int verbosity = 0);

  // Attribute accessors
  int get_size() { return size; }
  int get_no_iter() { return no_iter; }
  Eigen::MatrixXd get_A() { return A.cast<double>(); }
  Eigen::VectorXd get_chiA() { return chiA.cast<double>(); }
  Eigen::VectorXd get_muD() { return muD.cast<double>(); }
  Eigen::VectorXd get_inv() { return inv.cast<double>(); }
  Eigen::VectorXd get_chev() { return chev.cast<double>(); }
  Eigen::MatrixXd get_D() { return D.cast<double>(); }
  Eigen::MatrixXd get_N() { return N.cast<double>(); }

  // Compute methods
  template <typename precision>
  JCDec& compute_chiA(bool round = false);
  JCDec& compute_chiA(std::string precsion = "double", bool round = true);
  JCDec& compute_muD(bool resultant = false);
  JCDec& compute_inv();
  JCDec& compute_chev();
  JCDec& compute_D_poly();
  JCDec& compute_D_mat();
  JCDec& compute_D(bool mat = false);
  JCDec& compute_N();
  JCDec& compute();

  // Check methods
  bool is_trivial();
  double check_cayleyhamilton();
  double check_nillpotency();
  double check_commutativity();

  // Meta attribute accessors
  double get_timing(std::string step_name = "");
  Eigen::RowVectorXd get_timings(int steps = 7);

  double get_maxBitsize(std::string step_name = "");
  Eigen::RowVectorXd get_maxBitsizes(int steps = 7);
};

#endif  // JORDANCHEVALLEY_HPP
