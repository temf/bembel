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
#ifndef BEMBEL_SRC_LINEAROPERATOR_DIFFERENTIALFORMENUM_HPP_
#define BEMBEL_SRC_LINEAROPERATOR_DIFFERENTIALFORMENUM_HPP_

namespace Bembel {
/**
 * \ingroup LinearOperator
 * \brief Provides information about the discrete space required for the
 * discretisation of a specific operator.
 **/
struct DifferentialForm {
  enum { Continuous = 0, DivConforming = 1, Discontinuous = 2 };
};

/**
 *    \ingroup LinearOperator
 *    \brief struct containing specifications on DifferentialForms, i.e., the
 * function spaces.
 **/
template <unsigned int DifferentialForm, typename Scalar>
struct DifferentialFormTraits {
  // empty, only to be used in its specialization
};
/**
 * \brief Specification of the DifferentialFormTraits for the continuous case.
 */
template <typename Scalar>
struct DifferentialFormTraits<DifferentialForm::Continuous, Scalar> {
  enum { FunctionSpaceVectorDimension = 1, FunctionSpaceOutputDimension = 1 };

  typedef Scalar FunctionSpaceValue;
};
/**
 * \brief Specification of the DifferentialFormTraits for the div conforming
 * case.
 */
template <typename Scalar>
struct DifferentialFormTraits<DifferentialForm::DivConforming, Scalar> {
  enum { FunctionSpaceVectorDimension = 2, FunctionSpaceOutputDimension = 3 };

  typedef Eigen::Matrix<Scalar, 3, 1> FunctionSpaceValue;
};
/**
 * \brief Specification of the DifferentialFormTraits for the discontinuous
 * case.
 */
template <typename Scalar>
struct DifferentialFormTraits<DifferentialForm::Discontinuous, Scalar> {
  enum { FunctionSpaceVectorDimension = 1, FunctionSpaceOutputDimension = 1 };

  typedef Scalar FunctionSpaceValue;
};

/**
 * \deprecated Use DifferentialFormTraits instead.
 */
template <unsigned int DF>
constexpr int getFunctionSpaceVectorDimension() {
  return DifferentialFormTraits<DF, double>::FunctionSpaceVectorDimension;
}

/**
 * \deprecated Use DifferentialFormTraits instead.
 */
template <unsigned int DF>
constexpr int getFunctionSpaceOutputDimension() {
  return DifferentialFormTraits<DF, double>::FunctionSpaceOutputDimension;
}

}  // namespace Bembel

#endif  // BEMBEL_SRC_LINEAROPERATOR_DIFFERENTIALFORMENUM_HPP_
