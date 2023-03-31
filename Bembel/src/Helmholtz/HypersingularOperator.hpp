// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
//
#ifndef BEMBEL_LINEAROPERATOR_HELMHOLTZ_HELMHOLTZHYPERSINGULAROPERATOR_H_
#define BEMBEL_LINEAROPERATOR_HELMHOLTZ_HELMHOLTZHYPERSINGULAROPERATOR_H_

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
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p1, const SurfacePoint &p2,
      Eigen::Matrix<
          typename LinearOperatorTraits<HelmholtzHypersingularOperator>::Scalar,
          Eigen::Dynamic, Eigen::Dynamic> *intval) const {
    auto polynomial_degree = super_space.get_polynomial_degree();
    auto polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);

    // get evaluation points on unit square
    auto s = p1.segment<2>(0);
    auto t = p2.segment<2>(0);

    // get quadrature weights
    auto ws = p1(2);
    auto wt = p2(2);

    // get points on geometry and tangential derivatives
    auto x_f = p1.segment<3>(3);
    auto x_f_dx = p1.segment<3>(6);
    auto x_f_dy = p1.segment<3>(9);
    auto y_f = p2.segment<3>(3);
    auto y_f_dx = p2.segment<3>(6);
    auto y_f_dy = p2.segment<3>(9);

    // compute surface measures from tangential derivatives
    auto x_n = x_f_dx.cross(x_f_dy);
    auto y_n = y_f_dx.cross(y_f_dy);

    // compute h
    auto h = 1. / (1 << super_space.get_refinement_level());  // h = 1 ./ (2^M)

    // evaluate kernel
    auto kernel = evaluateKernel(x_f, y_f);

    // integrand without basis functions
    auto integrandScalar = -kernel * x_n.dot(y_n) * wavenumber2_ * ws * wt;
    auto integrandCurl = kernel * x_n.norm() * y_n.norm() * ws * wt / h / h;

    // multiply basis functions with integrand and add to intval, this is an
    // efficient implementation of
    super_space.addScaledBasisInteraction(intval, integrandScalar, s, t);
    super_space.addScaledSurfaceCurlInteraction(intval, integrandCurl, p1, p2);

    return;
  }

  Eigen::Matrix<std::complex<double>, 3, 3> evaluateFMMInterpolation_impl(
      const SurfacePoint &p1, const SurfacePoint &p2) const {
    // get evaluation points on unit square
    auto s = p1.segment<2>(0);
    auto t = p2.segment<2>(0);

    // get points on geometry and tangential derivatives
    auto x_f = p1.segment<3>(3);
    auto x_f_dx = p1.segment<3>(6);
    auto x_f_dy = p1.segment<3>(9);
    auto y_f = p2.segment<3>(3);
    auto y_f_dx = p2.segment<3>(6);
    auto y_f_dy = p2.segment<3>(9);

    // compute surface measures from tangential derivatives
    auto x_n = x_f_dx.cross(x_f_dy);
    auto y_n = y_f_dx.cross(y_f_dy);

    // evaluate kernel
    auto kernel = evaluateKernel(x_f, y_f);

    // interpolation
    Eigen::Matrix<std::complex<double>, 3, 3> intval;
    intval.setZero();
    intval(0, 0) = -kernel * wavenumber2_ * x_n.dot(y_n);
    intval(1, 1) = kernel * x_f_dy.dot(y_f_dy);
    intval(1, 2) = -kernel * x_f_dy.dot(y_f_dx);
    intval(2, 1) = -kernel * x_f_dx.dot(y_f_dy);
    intval(2, 2) = kernel * x_f_dx.dot(y_f_dx);

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
#endif
