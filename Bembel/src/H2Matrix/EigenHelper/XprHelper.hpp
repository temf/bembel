// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_H2MATRIX_EIGENHELPER_XPRHELPER_H__
#define __BEMBEL_H2MATRIX_EIGENHELPER_XPRHELPER_H_

/**
 * This file mimicks Core/util/XprHelper.h of Eigen
 */

namespace Eigen {

namespace internal {

// remember: storage type is not necessarily storage type, but more a template
// parameter which handles the called specializations.

template <typename Functor>
struct cwise_promote_storage_type<H2, Dense, Functor> {
  typedef H2 ret;
};
template <typename Functor>
struct cwise_promote_storage_type<Dense, H2, Functor> {
  typedef H2 ret;
};
template <typename Functor>
struct cwise_promote_storage_type<H2, Sparse, Functor> {
  typedef H2 ret;
};
template <typename Functor>
struct cwise_promote_storage_type<Sparse, H2, Functor> {
  typedef H2 ret;
};

}  // end namespace internal

}  // end namespace Eigen

#endif