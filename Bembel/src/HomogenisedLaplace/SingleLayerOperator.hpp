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
template <>
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
class HomogenisedLaplaceSingleLayerOperator
    : public LinearOperatorBase<HomogenisedLaplaceSingleLayerOperator> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  /**
   * \brief Constructs an object initialising the coefficients and the degree
   *  via the static variable precision.
   */
  HomogenisedLaplaceSingleLayerOperator() {
    this->deg = getDegree(HomogenisedLaplaceSingleLayerOperator::precision);
    this->cs =
        getCoefficients(HomogenisedLaplaceSingleLayerOperator::precision);
  }

  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p1, const SurfacePoint &p2,
      Eigen::Matrix<typename LinearOperatorTraits<
                        HomogenisedLaplaceSingleLayerOperator>::Scalar,
                    Eigen::Dynamic, Eigen::Dynamic> *intval) const {
    // integrand without basis functions
    auto integrand = evaluateKernel(p1.get_f(), p2.get_f()) *
                     p1.get_surface_measure() * p2.get_surface_measure() *
                     p1.get_w() * p2.get_w();

    // multiply basis functions with integrand and add to intval
    super_space.addScaledBasisInteraction(intval, integrand, p1.get_xi(),
                                          p2.get_xi());

    return;
  }

  Eigen::Matrix<double, 1, 1> evaluateFMMInterpolation_impl(
      const SurfacePoint &p1, const SurfacePoint &p2) const {
    // interpolation
    Eigen::Matrix<double, 1, 1> intval;
    intval(0) = evaluateKernel(p1.get_f(), p2.get_f()) *
                p1.get_surface_measure() * p2.get_surface_measure();
    return intval;
  }

  /**
   * \brief Fundamental solution of the Homogenised Laplace problem
   */
  double evaluateKernel(const Eigen::Vector3d &x,
                        const Eigen::Vector3d &y) const {
    return k_mod(x - y) +
           evaluate_solid_sphericals(x - y, this->cs, this->deg, false);
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
