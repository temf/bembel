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

#ifndef BEMBEL_SRC_LINEAROPERATOR_LOCALOPERATORBASE_HPP_
#define BEMBEL_SRC_LINEAROPERATOR_LOCALOPERATORBASE_HPP_

namespace Bembel {
/**
 *    \ingroup LocalOperator
 *    \brief local operator base class. this serves as a common interface for
 *           existing local operators
 **/
template <typename Derived>
struct LocalOperatorBase : public LinearOperatorBase<Derived> {
  // Constructors
  LocalOperatorBase() {}
  // the user has to provide the implementation of this function, which
  // is able to evaluate the integrand of the Galerkin formulation at a
  // quadrature point represented as a
  // Surface point [xi; h * w; Chi(xi); dsChi(xi); dtChi(xi)]
  using LinearOperatorBase<Derived>::evaluateIntegrand;
  // return the required quadrature degree for the far-field
  using LinearOperatorBase<Derived>::get_FarfieldQuadratureDegree;
  // pointer to the derived object
  using LinearOperatorBase<Derived>::derived;

 private:
  using LinearOperatorBase<Derived>::getNearfieldQuadratureDegree;
  using LinearOperatorBase<Derived>::evaluateFMMInterpolation;
};
}  // namespace Bembel
#endif  // BEMBEL_SRC_LINEAROPERATOR_LOCALOPERATORBASE_HPP_
