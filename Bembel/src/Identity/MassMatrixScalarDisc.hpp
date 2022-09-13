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

#ifndef BEMBEL_SRC_IDENTITY_MASSMATRIXSCALARDISC_HPP_
#define BEMBEL_SRC_IDENTITY_MASSMATRIXSCALARDISC_HPP_

namespace Bembel {
// forward declaration of class MassMatrixScalarDisc in order to define
// traits
class MassMatrixScalarDisc;

template <>
struct LinearOperatorTraits<MassMatrixScalarDisc> {
  typedef Eigen::VectorXd EigenType;
  typedef Eigen::VectorXd::Scalar Scalar;
  enum {
    OperatorOrder = 0,
    Form = DifferentialForm::Discontinuous,
    NumberOfFMMComponents = 0
  };
};

class MassMatrixScalarDisc : public IdentityOperatorBase<MassMatrixScalarDisc> {
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_IDENTITY_MASSMATRIXSCALARDISC_HPP_
