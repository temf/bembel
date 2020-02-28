// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#include <Bembel/AnsatzSpace>
#include "Test.hpp"

class TestOperatorDisc;
class TestOperatorDivC;

template <>
struct Bembel::LinearOperatorTraits<TestOperatorDisc> {
  typedef Eigen::VectorXd EigenType;
  typedef Eigen::VectorXd::Scalar Scalar;
  enum { OperatorOrder = 0, Form = DifferentialForm::Discontinuous };
};

template <>
struct Bembel::LinearOperatorTraits<TestOperatorDivC> {
  typedef Eigen::VectorXd EigenType;
  typedef Eigen::VectorXd::Scalar Scalar;
  enum { OperatorOrder = 0, Form = DifferentialForm::DivConforming };
};

int main() {
  using namespace Bembel;

  Test::TestGeometryWriter::writeScreen();

  Bembel::Geometry geometry("test_Screen.dat");
  assert(geometry.get_geometry().size() == 1);

  for (int p = 1; p < 5; ++p)
    for (int m = 0; m < 3; ++m)
      for (int k = 1; k < std::min(2, p + 1); ++k) {
        SuperSpace<TestOperatorDisc> superspace1(geometry, m, p);
        SuperSpace<TestOperatorDivC> superspace2(geometry, m, p);

        Projector<TestOperatorDisc> pro_sca(superspace1, k);
        Projector<TestOperatorDivC> pro_vec(superspace2, k);
        Eigen::SparseMatrix<double> mat_sca = pro_sca.get_projection_matrix();
        Eigen::SparseMatrix<double> mat_vec = pro_vec.get_projection_matrix();

        // These asserts check if the numbers of degrees of freedom are as
        // expected
        BEMBEL_TEST_IF(mat_sca.rows() ==
                       (p + 1) * (p + 1) * (1 << (m)) * (1 << (m)));
        BEMBEL_TEST_IF(mat_sca.rows() * 2 == mat_vec.rows());
        BEMBEL_TEST_IF(mat_sca.cols() ==
                       ((p + (1 << m) + ((1 << m) - 1) * (k - 1)) *
                        (p + (1 << m) + ((1 << m) - 1) * (k - 1))));
        BEMBEL_TEST_IF(mat_vec.cols() ==
                       2 * ((p + (1 << m) + ((1 << m) - 1) * (k - 1)) *
                            (p - 1 + (1 << m) + ((1 << m) - 1) * (k - 1))));

        Eigen::VectorXd coef_sca = Eigen::VectorXd::Ones(mat_sca.cols());
        Eigen::VectorXd coef_vec = Eigen::VectorXd::Ones(mat_vec.cols());

        Eigen::VectorXd big_coef_sca = mat_sca * coef_sca;
        Eigen::VectorXd big_coef_vec = mat_vec * coef_vec;

        // Both the superspace and the spline space offer a partition of unity
        // with 1-coefficients
        for (int i = 0; i < mat_sca.rows(); ++i)
          BEMBEL_TEST_IF((std::abs(big_coef_sca(i) - 1.0)) <
                         Test::Constants::coefficient_accuracy);

        // Both the superspace and the spline space offer a partition of unity
        // with 1-coefficients
        for (int i = 0; i < mat_vec.rows(); ++i)
          BEMBEL_TEST_IF((std::abs(big_coef_vec(i) - 1.0)) <
                         Test::Constants::coefficient_accuracy);

        // From here on we check the scalar projector against the deBoor
        // Formula. First, we fix some arbitrary coefficients.
        for (int i = 0; i < mat_sca.cols(); ++i)
          coef_sca(i) = (i % (p + 1 + m)) * 1.0;

        // This should default to the space that corresponds to the big space up
        // to ordering
        Projector<TestOperatorDisc> pro_big(superspace1, 1000);
        Eigen::SparseMatrix<double> mat_big = pro_big.get_projection_matrix();

        // Now, we blow up the coefficients.
        big_coef_sca = mat_big.transpose() * mat_sca * coef_sca;

        // Here, we construct the knot vector that should represent the small
        // space
        auto small_knot_vec =
            Bembel::Spl::MakeUniformKnotVector(p + 1, (1 << m) - 1, k);

        // Here, we construct the knot vector that should represent the
        // large space
        auto large_knot_vec =
            Bembel::Spl::MakeUniformKnotVector(p + 1, (1 << m) - 1, p + 1);

        const std::vector<double> test_points = {0.0, 0.1, 0.2, 0.3, 0.4,
                                                 0.5, 0.6, 0.7, 0.8, 1.0};

        // The deBoor Formula needs the coefficients to be in a matrix
        Eigen::MatrixXd coef_mat_large((p + 1) * (1 << (m)),
                                       (p + 1) * (1 << (m)));
        Eigen::MatrixXd coef_mat_small(
            (p + (1 << m) + ((1 << m) - 1) * (k - 1)),
            (p + (1 << m) + ((1 << m) - 1) * (k - 1)));

        for (int i = 0; i < coef_mat_large.rows() * coef_mat_large.cols(); ++i)
          coef_mat_large(i) = big_coef_sca(i);

        for (int i = 0; i < coef_mat_small.rows() * coef_mat_small.cols(); ++i)
          coef_mat_small(i) = coef_sca(i);

        // Evaluation via Recursion formula
        Eigen::MatrixXd result1 =
            Bembel::Spl::DeBoorTP(coef_mat_small, small_knot_vec,
                                  small_knot_vec, test_points, test_points);
        Eigen::MatrixXd result2 =
            Bembel::Spl::DeBoorTP(coef_mat_large, large_knot_vec,
                                  large_knot_vec, test_points, test_points);
        // Testing of the result
        BEMBEL_TEST_IF((result1 - result2).norm() <
                       Test::Constants::coefficient_accuracy);
      }

  return 0;
}
