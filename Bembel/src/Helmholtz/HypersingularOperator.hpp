// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2024 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
//
#ifndef BEMBEL_SRC_HELMHOLTZ_HYPERSINGULAROPERATOR_HPP_
#define BEMBEL_SRC_HELMHOLTZ_HYPERSINGULAROPERATOR_HPP_

namespace Bembel {
// forward declaration of class HelmholtzHypersingularOperator in order to
// define traits
class HelmholtzHypersingularOperator;

template <>
struct LinearOperatorTraits<HelmholtzHypersingularOperator> {
  typedef Eigen::VectorXcd EigenType;
  typedef Eigen::VectorXcd::Scalar Scalar;
  enum {
    OperatorOrder = 1,
    Form = DifferentialForm::Continuous,
    NumberOfFMMComponents = 3
  };
};

/**
 * \ingroup Helmholtz
 */
class HelmholtzHypersingularOperator
    : public LinearOperatorBase<HelmholtzHypersingularOperator> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  HelmholtzHypersingularOperator() {}
  template <class T>
  void evaluateIntegrand_impl(const T &super_space, const SurfacePoint &p1,
                              const SurfacePoint &p2,
                              Eigen::MatrixXcd *intval) const {
    // compute surface measures from tangential derivatives
    Eigen::Vector3d x_n = p1.get_normal();
    Eigen::Vector3d y_n = p2.get_normal();

    // compute h
    double h =
        1. / (1 << super_space.get_refinement_level());  // h = 1 ./ (2^M)

    // evaluate kernel
    std::complex<double> kernel = evaluateKernel(p1.get_f(), p2.get_f());

    // integrand without basis functions
    std::complex<double> integrandScalar =
        -kernel * x_n.dot(y_n) * wavenumber2_ * p1.get_w() * p2.get_w();
    std::complex<double> integrandCurl =
        kernel * x_n.norm() * y_n.norm() * p1.get_w() * p2.get_w() / h / h;

    // multiply basis functions with integrand and add to intval, this is an
    // efficient implementation of
    super_space.addScaledBasisInteraction(intval, integrandScalar, p1.get_xi(),
                                          p2.get_xi());
    super_space.addScaledSurfaceCurlInteraction(intval, integrandCurl, p1, p2);

    return;
  }

  Eigen::Matrix<std::complex<double>, 3, 3> evaluateFMMInterpolation_impl(
      const SurfacePoint &p1, const SurfacePoint &p2) const {
    // evaluate kernel
    std::complex<double> kernel = evaluateKernel(p1.get_f(), p2.get_f());

    // interpolation
    Eigen::Matrix<std::complex<double>, 3, 3> intval;
    intval.setZero();
    intval(0, 0) =
        -kernel * wavenumber2_ * p1.get_normal().dot(p2.get_normal());
    intval(1, 1) = kernel * p1.get_f_dy().dot(p2.get_f_dy());
    intval(1, 2) = -kernel * p1.get_f_dy().dot(p2.get_f_dx());
    intval(2, 1) = -kernel * p1.get_f_dx().dot(p2.get_f_dy());
    intval(2, 2) = kernel * p1.get_f_dx().dot(p2.get_f_dx());

    return intval;
  }

  /**
   * \brief Fundamental solution of Helmholtz problem
   */
  std::complex<double> evaluateKernel(const Eigen::Vector3d &x,
                                      const Eigen::Vector3d &y) const {
    auto r = (x - y).norm();
    return std::exp(-std::complex<double>(0., 1.) * wavenumber_ * r) / 4. /
           BEMBEL_PI / r;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    setters
  //////////////////////////////////////////////////////////////////////////////
  void set_wavenumber(std::complex<double> wavenumber) {
    wavenumber_ = wavenumber;
    wavenumber2_ = wavenumber_ * wavenumber_;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    getters
  //////////////////////////////////////////////////////////////////////////////
  std::complex<double> get_wavenumber() { return wavenumber_; }

 private:
  std::complex<double> wavenumber_;
  std::complex<double> wavenumber2_;
};

/**
 * \brief The hypersingular operator requires a special treatment of the
 * moment matrices of the FMM due to the involved derivatives on the ansatz
 * functions.
 */
template <typename InterpolationPoints>
struct H2Multipole::Moment2D<InterpolationPoints,
                             HelmholtzHypersingularOperator> {
  static std::vector<Eigen::MatrixXd> compute2DMoment(
      const SuperSpace<HelmholtzHypersingularOperator> &super_space,
      const int cluster_level, const int cluster_refinements,
      const int number_of_points) {
    Eigen::MatrixXd moment = moment2DComputer<
        Moment1D<InterpolationPoints, HelmholtzHypersingularOperator>,
        Moment1D<InterpolationPoints, HelmholtzHypersingularOperator>>(
        super_space, cluster_level, cluster_refinements, number_of_points);
    Eigen::MatrixXd moment_dx = moment2DComputer<
        Moment1DDerivative<InterpolationPoints, HelmholtzHypersingularOperator>,
        Moment1D<InterpolationPoints, HelmholtzHypersingularOperator>>(
        super_space, cluster_level, cluster_refinements, number_of_points);
    Eigen::MatrixXd moment_dy = moment2DComputer<
        Moment1D<InterpolationPoints, HelmholtzHypersingularOperator>,
        Moment1DDerivative<InterpolationPoints,
                           HelmholtzHypersingularOperator>>(
        super_space, cluster_level, cluster_refinements, number_of_points);

    Eigen::MatrixXd moment_total(
        moment.rows() + moment_dx.rows() + moment_dy.rows(), moment_dx.cols());
    moment_total << moment, moment_dx, moment_dy;

    std::vector<Eigen::MatrixXd> vector_of_moments;
    vector_of_moments.push_back(moment_total);

    return vector_of_moments;
  }
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_HELMHOLTZ_HYPERSINGULAROPERATOR_HPP_
