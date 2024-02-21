// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2022 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_SRC_LINEAROPERATOR_DUMMY_DUMMYOPERATOR_HPP_
#define BEMBEL_SRC_LINEAROPERATOR_DUMMY_DUMMYOPERATOR_HPP_

namespace Bembel {

class DummyOperator;

template <>
struct LinearOperatorTraits<DummyOperator> {
  typedef Eigen::VectorXd EigenType;
  typedef Eigen::VectorXd::Scalar Scalar;
  enum {
    OperatorOrder = 0,
    Form = DifferentialForm::Discontinuous,
    NumberOfFMMComponents = 1
  };
};

// forward declaration of class DummyOperator in order to define traits
// define some default test functiom
std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)>
    DummyOperator_test_function =
        [](const Eigen::Vector2d &x, const Eigen::Vector2d &y) { return 1.; };

/**
 * \ingroup DummyOperator
 * \brief This class provides a dummy specialization of the LinearOperator and
 * corresponding Traits for testing and debugging
 */
class DummyOperator : public LinearOperatorBase<DummyOperator> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  /**
   * \brief Default constructor
   *
   * The default test function is a constant 1.
   */
  DummyOperator() { test_func_ = DummyOperator_test_function; }
  /**
   * \brief Constructor with a given test function.
   */
  DummyOperator(
      std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)>
          test_func) {
    test_func_ = test_func;
  }
  /**
   * \brief Implements the integration routine.
   *
   * \param super_space SuperSpace specified with template parameter.
   * \param p1 SurfacePoint for evaluating the first integral.
   * \param p2 SurfacePoint for evaluating the second integral.
   * \param intval Matrix for the computed integral.
   */
  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p1, const SurfacePoint &p2,
      Eigen::Matrix<typename LinearOperatorTraits<DummyOperator>::Scalar,
                    Eigen::Dynamic, Eigen::Dynamic> *intval) const {
    (*intval)(0, 0) +=
        test_func_(p1.segment(3, 2), p2.segment(3, 2)) * p1(2) * p2(2);
    return;
  }

 private:
  std::function<double(const Eigen::Vector2d &, const Eigen::Vector2d &)>
      test_func_;
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_LINEAROPERATOR_DUMMY_DUMMYOPERATOR_HPP_
