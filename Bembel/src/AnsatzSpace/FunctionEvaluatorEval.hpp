// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_ANSATZSPACE_FUNCTIONEVALUATOREVAL_H_
#define BEMBEL_ANSATZSPACE_FUNCTIONEVALUATOREVAL_H_

namespace Bembel {

template <typename Scalar, unsigned int DF, typename LinOp>
struct FunctionEvaluatorEval {};

// continuous
template <typename Scalar, typename LinOp>
struct FunctionEvaluatorEval<Scalar, DifferentialForm::Continuous, LinOp> {
  Eigen::Matrix<Scalar,
                getFunctionSpaceOutputDimension<DifferentialForm::Continuous>(),
                1>
  eval(const SuperSpace<LinOp> &super_space,
       const int polynomial_degree_plus_one_squared,
       const ElementTreeNode &element, const SurfacePoint &p,
       const Eigen::Matrix<
           Scalar, Eigen::Dynamic,
           getFunctionSpaceVectorDimension<DifferentialForm::Continuous>()>
           &coeff) const {
    auto s = p.segment<2>(0);
    return coeff.transpose() * super_space.basis(s) / element.get_h();
  };
};

// div-conforming
template <typename Scalar, typename LinOp>
struct FunctionEvaluatorEval<Scalar, DifferentialForm::DivConforming, LinOp> {
  Eigen::Matrix<
      Scalar,
      getFunctionSpaceOutputDimension<DifferentialForm::DivConforming>(), 1>
  eval(const SuperSpace<LinOp> &super_space,
       const int polynomial_degree_plus_one_squared,
       const ElementTreeNode &element, const SurfacePoint &p,
       const Eigen::Matrix<
           Scalar, Eigen::Dynamic,
           getFunctionSpaceVectorDimension<DifferentialForm::DivConforming>()>
           &coeff) const {
    auto s = p.segment<2>(0);
    auto h = element.get_h();
    auto x_f_dx = p.segment<3>(6);
    auto x_f_dy = p.segment<3>(9);
    Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar, Eigen::Dynamic,
                  1>
        tangential_coefficients = coeff.transpose() * super_space.basis(s);
    return (x_f_dx * tangential_coefficients(0) +
            x_f_dy * tangential_coefficients(1)) /
           h;
  };

  Scalar evalDiv(
      const SuperSpace<LinOp> &super_space,
      const int polynomial_degree_plus_one_squared,
      const ElementTreeNode &element, const SurfacePoint &p,
      const Eigen::Matrix<
          Scalar, Eigen::Dynamic,
          getFunctionSpaceVectorDimension<DifferentialForm::DivConforming>()>
          &coeff) const {
    auto s = p.segment<2>(0);
    auto h = element.get_h();
    Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar, Eigen::Dynamic,
                  1>
        phiPhiVec_dx = super_space.basisDx(s);
    Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar, Eigen::Dynamic,
                  1>
        phiPhiVec_dy = super_space.basisDy(s);
    return (phiPhiVec_dx.dot(coeff.col(0)) + phiPhiVec_dy.dot(coeff.col(1))) /
           h / h;
  };
};

// discontinuous
template <typename Scalar, typename LinOp>
struct FunctionEvaluatorEval<Scalar, DifferentialForm::Discontinuous, LinOp> {
  Eigen::Matrix<
      Scalar,
      getFunctionSpaceOutputDimension<DifferentialForm::Discontinuous>(), 1>
  eval(const SuperSpace<LinOp> &super_space,
       const int polynomial_degree_plus_one_squared,
       const ElementTreeNode &element, const SurfacePoint &p,
       const Eigen::Matrix<
           Scalar, Eigen::Dynamic,
           getFunctionSpaceVectorDimension<DifferentialForm::Discontinuous>()>
           &coeff) const {
    auto s = p.segment<2>(0);
    return coeff.transpose() * super_space.basis(s) / element.get_h();
  };
};
}  // namespace Bembel
#endif
