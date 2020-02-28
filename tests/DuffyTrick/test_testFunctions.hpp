// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
//
#ifndef BEMBEL_TEST_DUFFYTRICK_TESTFUNCTIONS_H_
#define BEMBEL_TEST_DUFFYTRICK_TESTFUNCTIONS_H_

namespace Test {
namespace DuffyTrick {

static const double sigma = 1e-3;

/**
 *  \brief dependent on sigma, this is an almost singular function with an
 *         antiderivative which can easily be computed
 **/
double testFunction1_2D(const Eigen::Vector2d &x) {
  return 1. / (sigma + BEMBEL_SQUARED_(x(0) - x(1)));
}

double testFuncion1_2DAntiDer(const Eigen::Vector2d &x) {
  return 0.5 * std::log(sigma + BEMBEL_SQUARED_(x(0) - x(1))) +
         (x(1) - x(0)) / sqrt(sigma) * std::atan((x(0) - x(1)) / sqrt(sigma));
}
////////////////////////////////////////////////////////////////////////////////
/**
 *  \brief analytic integral of the the 2D testfunction wrt a rectangle
 **/
double testFunction2DIntegral(const Eigen::Matrix2d &axis) {
  double retval = 0;
  Eigen::Vector2d curr_pt;
  Eigen::Vector4d intvals;
  curr_pt << axis(1, 0), axis(1, 1);
  retval += testFuncion1_2DAntiDer(curr_pt);
  intvals(0) = testFuncion1_2DAntiDer(curr_pt);
  curr_pt << axis(1, 0), axis(0, 1);
  retval -= testFuncion1_2DAntiDer(curr_pt);
  intvals(1) = testFuncion1_2DAntiDer(curr_pt);
  curr_pt << axis(0, 0), axis(1, 1);
  retval -= testFuncion1_2DAntiDer(curr_pt);
  intvals(2) = testFuncion1_2DAntiDer(curr_pt);
  curr_pt << axis(0, 0), axis(0, 1);
  retval += testFuncion1_2DAntiDer(curr_pt);
  intvals(3) = testFuncion1_2DAntiDer(curr_pt);
  return retval;
}
////////////////////////////////////////////////////////////////////////////////
// we use the coupling function in 2D to generate a test function
// in 4D which couples to elements
double testFunction4D(const Eigen::Vector2d &x, const Eigen::Vector2d &y) {
  Eigen::Vector2d coupling;
  coupling << x(1), y(0);
  return std::exp(-x(0) - 2 * y(1)) * testFunction1_2D(coupling);
}
// analytic integral of the 4D function
double testFunction4DIntegral(const Eigen::MatrixXd &axis) {
  return (std::exp(-axis(0, 0)) - std::exp(-axis(1, 0))) *
         testFunction2DIntegral(axis.block(0, 1, 2, 2)) *
         (0.5 * std::exp(-2 * axis(0, 3)) - 0.5 * std::exp(-2 * axis(1, 3)));
}
////////////////////////////////////////////////////////////////////////////////
/**
 *  \brief this function uses a simple tensor product quadrature, we test
 *  it with an analytic function
 **/
std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)>
    integrate0_test_function =
        [](const Eigen::Vector2d &x, const Eigen::Vector2d &y) {
          return std::sin(BEMBEL_PI * x(0)) * std::sin(2 * BEMBEL_PI * y(0)) *
                 std::exp(-x(1)) * std::exp(y(1));
        };
std::function<double(const Eigen::MatrixXd &)>
    integrate0_test_function_integral = [](const Eigen::MatrixXd &axis) {
      return (std::cos(BEMBEL_PI * axis(0, 0)) -
              std::cos(BEMBEL_PI * axis(1, 0))) /
             BEMBEL_PI * 0.5 *
             (std::cos(2 * BEMBEL_PI * axis(0, 2)) -
              std::cos(2 * BEMBEL_PI * axis(1, 2))) /
             BEMBEL_PI * (std::exp(-axis(0, 1)) - std::exp(-axis(1, 1))) *
             (std::exp(axis(1, 3)) - std::exp(axis(0, 3)));
    };
/**
 *  \brief this function uses a simple tensor product quadrature, we test
 *  it with an analytic function
 **/
std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)>
    integrate1_test_function =
        [](const Eigen::Vector2d &x, const Eigen::Vector2d &y) {
          return std::sin(BEMBEL_PI * x(0)) * std::sin(2 * BEMBEL_PI * y(0)) *
                 std::exp(-x(1)) * std::exp(y(1));
        };
std::function<double(const Eigen::MatrixXd &)>
    integrate1_test_function_integral = [](const Eigen::MatrixXd &axis) {
      return (std::cos(BEMBEL_PI * axis(0, 0)) -
              std::cos(BEMBEL_PI * axis(1, 0))) /
             BEMBEL_PI * 0.5 *
             (std::cos(2 * BEMBEL_PI * axis(0, 2)) -
              std::cos(2 * BEMBEL_PI * axis(1, 2))) /
             BEMBEL_PI * (std::exp(-axis(0, 1)) - std::exp(-axis(1, 1))) *
             (std::exp(axis(1, 3)) - std::exp(axis(0, 3)));
    };

std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)>
    integrate2_test_function =
        [](const Eigen::Vector2d &x, const Eigen::Vector2d &y) {
#ifdef BEMBEL_TEST_DUFFYTRICK_USE_ANALYTIC_FUNCTION_
          return std::sin(BEMBEL_PI * x(0)) * std::sin(2 * BEMBEL_PI * y(0)) *
                 std::exp(-x(1)) * std::exp(y(1));
#else
          return testFunction4D(x, y);
#endif
        };
std::function<double(const Eigen::MatrixXd &)>
    integrate2_test_function_integral = [](const Eigen::MatrixXd &axis) {
#ifdef BEMBEL_TEST_DUFFYTRICK_USE_ANALYTIC_FUNCTION_
      return (std::cos(BEMBEL_PI * axis(0, 0)) -
              std::cos(BEMBEL_PI * axis(1, 0))) /
             BEMBEL_PI * 0.5 *
             (std::cos(2 * BEMBEL_PI * axis(0, 2)) -
              std::cos(2 * BEMBEL_PI * axis(1, 2))) /
             BEMBEL_PI * (std::exp(-axis(0, 1)) - std::exp(-axis(1, 1))) *
             (std::exp(axis(1, 3)) - std::exp(axis(0, 3)));
#else
      return testFunction4DIntegral(axis);
#endif
    };

std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)>
    integrate3_test_function =
        [](const Eigen::Vector2d &x, const Eigen::Vector2d &y) {
#ifdef BEMBEL_TEST_DUFFYTRICK_USE_ANALYTIC_FUNCTION_
          return std::sin(BEMBEL_PI * x(0)) * std::sin(2 * BEMBEL_PI * y(0)) *
                 std::exp(-x(1)) * std::exp(y(1));
#else
          return testFunction4D(x, y);
#endif
        };
std::function<double(const Eigen::MatrixXd &)>
    integrate3_test_function_integral = [](const Eigen::MatrixXd &axis) {
#ifdef BEMBEL_TEST_DUFFYTRICK_USE_ANALYTIC_FUNCTION_
      return (std::cos(BEMBEL_PI * axis(0, 0)) -
              std::cos(BEMBEL_PI * axis(1, 0))) /
             BEMBEL_PI * 0.5 *
             (std::cos(2 * BEMBEL_PI * axis(0, 2)) -
              std::cos(2 * BEMBEL_PI * axis(1, 2))) /
             BEMBEL_PI * (std::exp(-axis(0, 1)) - std::exp(-axis(1, 1))) *
             (std::exp(axis(1, 3)) - std::exp(axis(0, 3)));
#else
      return testFunction4DIntegral(axis);
#endif
    };

std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)>
    integrate4_test_function =
        [](const Eigen::Vector2d &x, const Eigen::Vector2d &y) {
#ifdef BEMBEL_TEST_DUFFYTRICK_USE_ANALYTIC_FUNCTION_
          return std::sin(BEMBEL_PI * x(0)) * std::sin(2 * BEMBEL_PI * y(0)) *
                 std::exp(-x(1)) * std::exp(y(1));
#else
          return testFunction4D(x, y);
#endif
        };
std::function<double(const Eigen::MatrixXd &)>
    integrate4_test_function_integral = [](const Eigen::MatrixXd &axis) {
#ifdef BEMBEL_TEST_DUFFYTRICK_USE_ANALYTIC_FUNCTION_
      return (std::cos(BEMBEL_PI * axis(0, 0)) -
              std::cos(BEMBEL_PI * axis(1, 0))) /
             BEMBEL_PI * 0.5 *
             (std::cos(2 * BEMBEL_PI * axis(0, 2)) -
              std::cos(2 * BEMBEL_PI * axis(1, 2))) /
             BEMBEL_PI * (std::exp(-axis(0, 1)) - std::exp(-axis(1, 1))) *
             (std::exp(axis(1, 3)) - std::exp(axis(0, 3)));
#else
      return testFunction4DIntegral(axis);
#endif
    };
}  // namespace DuffyTrick
}  // namespace Test
#endif
