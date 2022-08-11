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
#ifndef BEMBEL_SRC_POTENTIAL_POTENTIAL_HPP_
#define BEMBEL_SRC_POTENTIAL_POTENTIAL_HPP_

#include <Eigen/Dense>

#include "../LinearOperator/DifferentialFormEnum.hpp"
#include "../LinearOperator/LinearOperatorTraits.hpp"
#include "../util/Macros.hpp"

namespace Bembel {
/**
 *    \ingroup Potential
 *    \brief struct containing specifications on the functional
 *           has to be specialized or derived for any particular operator
 *           under consideration
 **/
template <typename Derived>
struct PotentialTraits {
  enum { YOU_DID_NOT_SPECIFY_POTENTIAL_TRAITS = 1 };
};

template <typename S, typename T>
struct PotentialReturnScalar {
  enum { RETURN_TYPE_ONLY_SPECIFIED_FOR_DOUBLE_OR_COMPLEX_DOUBLE = 1 };
};
template <>
struct PotentialReturnScalar<double, double> {
  typedef double Scalar;
};
template <>
struct PotentialReturnScalar<std::complex<double>, double> {
  typedef std::complex<double> Scalar;
};
template <>
struct PotentialReturnScalar<double, std::complex<double>> {
  typedef std::complex<double> Scalar;
};
template <>
struct PotentialReturnScalar<std::complex<double>, std::complex<double>> {
  typedef std::complex<double> Scalar;
};

/**
 *    \brief functional base class. this serves as a common interface for
 *           existing functionals
 **/
template <typename Derived, typename LinOp>
struct PotentialBase {
  // Constructors
  PotentialBase() {}

  // the user has to provide the implementation of this function, which
  // tells
  // is able to evaluate the integrand of the Galerkin formulation in a
  // pair
  // of quadrature points represented as a
  // Surface point [xi; w; Chi(xi); dsChi(xi); dtChi(xi)]
  template <typename T>
  Eigen::Matrix<typename PotentialReturnScalar<
                    typename LinearOperatorTraits<LinOp>::Scalar,
                    typename PotentialTraits<Derived>::Scalar>::Scalar,
                1, PotentialTraits<Derived>::OutputSpaceDimension>
  evaluateIntegrand(
      const AnsatzSpace<LinOp> &ansatz_space,
      const Eigen::Matrix<
          typename LinearOperatorTraits<LinOp>::Scalar, Eigen::Dynamic,
          getFunctionSpaceVectorDimension<LinearOperatorTraits<LinOp>::Form>()>
          &coeff,
      const Eigen::VectorXd &point, const SurfacePoint &p) const {
    return static_cast<const Derived *>(this)->evaluatePotential_impl(
        ansatz_space, coeff, point, p);
  }
  // pointer to the derived object
  Derived &derived() { return *static_cast<Derived *>(this); }
  // const pointer to the derived object
  const Derived &derived() const { return *static_cast<const Derived *>(this); }
};
}  // namespace Bembel
#endif  // BEMBEL_SRC_POTENTIAL_POTENTIAL_HPP_
