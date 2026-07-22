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

#ifndef BEMBEL_SRC_IDENTITY_IDENTITYOPERATORBASE_HPP_
#define BEMBEL_SRC_IDENTITY_IDENTITYOPERATORBASE_HPP_

namespace Bembel {

/**
 * \ingroup LocalOperator
 * \brief This class is the base for all mass matrices.
 */
template <typename Derived>
class IdentityOperatorBase : public LocalOperatorBase<Derived> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  IdentityOperatorBase() {}
  template <class T>
  void evaluateIntegrand_impl(const T &super_space, const SurfacePoint &p1,
                              const SurfacePoint &p2,
                              Eigen::MatrixXd *intval) const {
    // integrand without basis functions
    const auto integrand = p1.get_surface_measure() * p1.get_w();
    super_space.addScaledBasisInteraction(intval, integrand, p1.get_xi(),
                                          p1.get_xi());
    return;
  }
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_IDENTITY_IDENTITYOPERATORBASE_HPP_
