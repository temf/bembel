// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_LINEARFORM_LINEARFORM_H_
#define BEMBEL_LINEARFORM_LINEARFORM_H_

namespace Bembel {
/**
 *    \ingroup LinearForm
 *    \brief This class needs to be specialized, such that key traits for user
 *defined LinearForms are available.
 **/
template <typename Derived>
struct LinearFormTraits {
  enum { YOU_DID_NOT_SPECIFY_LINEARFORM_TRAITS = 1 };
};

/**
 *  \ingroup LinearForm
 *  \brief This class provides a blueprint for the class that needs to be
 * specialized for assembly of the right hand side of the linear system.
 */
template <typename Derived, typename Scalar>
struct LinearFormBase {
  // Constructors
  LinearFormBase(){};

  // the user has to provide the implementation of this function, which
  // tells
  // is able to evaluate the integrand of the Galerkin formulation in a
  // pair
  // of quadrature points represented as a
  // Surface point [xi; w; Chi(xi); dsChi(xi); dtChi(xi)]
  template <class T>
  void evaluateIntegrand(
      const T &super_space, const SurfacePoint &p,
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *intval) const {
    static_cast<const Derived *>(this)->evaluateLinearForm_impl(super_space, p,
                                                                intval);
    return;
  }
  // pointer to the derived object
  Derived &derived() { return *static_cast<Derived *>(this); }
  // const pointer to the derived object
  const Derived &derived() const { return *static_cast<const Derived *>(this); }
};
}  // namespace Bembel
#endif
