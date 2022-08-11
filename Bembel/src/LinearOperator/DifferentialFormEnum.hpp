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
 *
 **/
struct DifferentialForm {
  enum { Continuous = 0, DivConforming = 1, Discontinuous = 2 };
};

template <unsigned int DF>
constexpr int getFunctionSpaceVectorDimension() {
  return -1;
}

template <unsigned int DF>
constexpr int getFunctionSpaceOutputDimension() {
  return -1;
}

template <>
constexpr int getFunctionSpaceVectorDimension<DifferentialForm::Continuous>() {
  return 1;
}

template <>
constexpr int getFunctionSpaceOutputDimension<DifferentialForm::Continuous>() {
  return 1;
}

template <>
constexpr int
getFunctionSpaceVectorDimension<DifferentialForm::DivConforming>() {
  return 2;
}

template <>
constexpr int
getFunctionSpaceOutputDimension<DifferentialForm::DivConforming>() {
  return 3;
}

template <>
constexpr int
getFunctionSpaceVectorDimension<DifferentialForm::Discontinuous>() {
  return 1;
}

template <>
constexpr int
getFunctionSpaceOutputDimension<DifferentialForm::Discontinuous>() {
  return 1;
}

}  // namespace Bembel

#endif  // BEMBEL_SRC_LINEAROPERATOR_DIFFERENTIALFORMENUM_HPP_
