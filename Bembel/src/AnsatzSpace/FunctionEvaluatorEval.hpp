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
#ifndef BEMBEL_SRC_ANSATZSPACE_FUNCTIONEVALUATOREVAL_HPP_
#define BEMBEL_SRC_ANSATZSPACE_FUNCTIONEVALUATOREVAL_HPP_

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
    return coeff.transpose() * super_space.basis(p.get_xi()) / element.get_h();
  }
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
    double h = element.get_h();
    Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar, Eigen::Dynamic,
                  1>
        tangential_coefficients =
            coeff.transpose() * super_space.basis(p.get_xi());
    return (p.get_f_dx() * tangential_coefficients(0) +
            p.get_f_dy() * tangential_coefficients(1)) /
           h;
  }

  Scalar evalDiv(
      const SuperSpace<LinOp> &super_space,
      const int polynomial_degree_plus_one_squared,
      const ElementTreeNode &element, const SurfacePoint &p,
      const Eigen::Matrix<
          Scalar, Eigen::Dynamic,
          getFunctionSpaceVectorDimension<DifferentialForm::DivConforming>()>
          &coeff) const {
    double h = element.get_h();
    Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar, Eigen::Dynamic,
                  1>
        phiPhiVec_dx = super_space.basisDx(p.get_xi());
    Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar, Eigen::Dynamic,
                  1>
        phiPhiVec_dy = super_space.basisDy(p.get_xi());
    return (phiPhiVec_dx.dot(coeff.col(0)) + phiPhiVec_dy.dot(coeff.col(1))) /
           h / h;
  }
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
    return coeff.transpose() * super_space.basis(p.get_xi()) / element.get_h();
  }
};
}  // namespace Bembel
#endif  // BEMBEL_SRC_ANSATZSPACE_FUNCTIONEVALUATOREVAL_HPP_
