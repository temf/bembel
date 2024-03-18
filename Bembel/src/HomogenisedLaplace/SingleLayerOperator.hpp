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

#ifndef BEMBEL_SRC_HOMOGENISEDLAPLACE_SINGLELAYEROPERATOR_HPP_
#define BEMBEL_SRC_HOMOGENISEDLAPLACE_SINGLELAYEROPERATOR_HPP_

namespace Bembel {
// forward declaration of class HomogenisedLaplaceSingleLayerOperator
// in order to define traits
class HomogenisedLaplaceSingleLayerOperator;

/**
 * \brief Specification of the LinerOperatorTraits for the Homogenised Laplace.
 */
template<>
struct LinearOperatorTraits<HomogenisedLaplaceSingleLayerOperator> {
  typedef Eigen::VectorXd EigenType;
  typedef Eigen::VectorXd::Scalar Scalar;
  enum {
    OperatorOrder = -1,
    Form = DifferentialForm::Discontinuous,
    NumberOfFMMComponents = 1
  };
};

/**
 * \ingroup HomogenisedLaplace
 * \brief This class implements the specification of the integration for the
 * single layer operator for the homogenised Laplace.
 */
class HomogenisedLaplaceSingleLayerOperator : public LinearOperatorBase<
    HomogenisedLaplaceSingleLayerOperator> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  /**
   * \brief Constructs an object initialising the coefficients and the degree
   *  via the static variable precision.
   */
  HomogenisedLaplaceSingleLayerOperator() {
    this->deg = getDegree(HomogenisedLaplaceSingleLayerOperator::precision);
    this->cs = getCoefficients(
        HomogenisedLaplaceSingleLayerOperator::precision);
  }

  template<class T>
  void evaluateIntegrand_impl(const T &super_space, const SurfacePoint &p1,
      const SurfacePoint &p2,
      Eigen::Matrix<
          typename LinearOperatorTraits<HomogenisedLaplaceSingleLayerOperator
          >::Scalar, Eigen::Dynamic, Eigen::Dynamic> *intval) const {
    auto polynomial_degree = super_space.get_polynomial_degree();
    auto polynomial_degree_plus_one_squared = (polynomial_degree + 1)
        * (polynomial_degree + 1);

    // get evaluation points on unit square
    auto s = p1.segment < 2 > (0);
    auto t = p2.segment < 2 > (0);

    // get quadrature weights
    auto ws = p1(2);
    auto wt = p2(2);

    // get points on geometry and tangential derivatives
    auto x_f = p1.segment < 3 > (3);
    auto x_f_dx = p1.segment < 3 > (6);
    auto x_f_dy = p1.segment < 3 > (9);
    auto y_f = p2.segment < 3 > (3);
    auto y_f_dx = p2.segment < 3 > (6);
    auto y_f_dy = p2.segment < 3 > (9);

    // compute surface measures from tangential derivatives
    auto x_kappa = x_f_dx.cross(x_f_dy).norm();
    auto y_kappa = y_f_dx.cross(y_f_dy).norm();

    // integrand without basis functions
    auto integrand = evaluateKernel(x_f, y_f) * x_kappa * y_kappa * ws * wt;

    // multiply basis functions with integrand and add to intval, this is an
    // efficient implementation of
    // (*intval) += super_space.basisInteraction(s, t)
    //    * evaluateKernel(x_f, y_f) * x_kappa * y_kappa * ws * wt;
    super_space.addScaledBasisInteraction(intval, integrand, s, t);

    return;
  }

  Eigen::Matrix<double, 1, 1> evaluateFMMInterpolation_impl(
      const SurfacePoint &p1, const SurfacePoint &p2) const {
    // get evaluation points on unit square
    auto s = p1.segment < 2 > (0);
    auto t = p2.segment < 2 > (0);

    // get points on geometry and tangential derivatives
    auto x_f = p1.segment < 3 > (3);
    auto x_f_dx = p1.segment < 3 > (6);
    auto x_f_dy = p1.segment < 3 > (9);
    auto y_f = p2.segment < 3 > (3);
    auto y_f_dx = p2.segment < 3 > (6);
    auto y_f_dy = p2.segment < 3 > (9);

    // compute surface measures from tangential derivatives
    auto x_kappa = x_f_dx.cross(x_f_dy).norm();
    auto y_kappa = y_f_dx.cross(y_f_dy).norm();

    // interpolation
    Eigen::Matrix<double, 1, 1> intval;
    intval(0) = evaluateKernel(x_f, y_f) * x_kappa * y_kappa;

    return intval;
  }

  /**
   * \brief Fundamental solution of the Homogenised Laplace problem
   */
  double evaluateKernel(const Eigen::Vector3d &x,
      const Eigen::Vector3d &y) const {
    return k_mod(x - y)
        + evaluate_solid_sphericals(x - y, this->cs, this->deg, false);
  }

  /**
   * \brief sets the precision of the periodicity of the kernel
   */
  static void setPrecision(double p) {
    HomogenisedLaplaceSingleLayerOperator::precision = p;
  }

  /**
   * \brief returns the precision of the periodicity of the kernel
   */
  static double getPrecision() {
    return HomogenisedLaplaceSingleLayerOperator::precision;
  }

 private:
  /** The degree of the spherical harmonics expansion */
  unsigned int deg;
  /** The coefficients of the spherical harmonics expansion */
  Eigen::VectorXd cs;
  /** The precision of the periodicity of the kernel */
  static double precision;
};

double HomogenisedLaplaceSingleLayerOperator::precision = 0;

}  // namespace Bembel

#endif  // BEMBEL_SRC_HOMOGENISEDLAPLACE_SINGLELAYEROPERATOR_HPP_
