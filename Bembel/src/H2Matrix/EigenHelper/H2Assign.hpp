// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2ASSIGN_H__
#define __BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2ASSIGN_H__

// The contents of this file are a modification of the file
// SparseAssign.h from the Eigen library
namespace Eigen {

namespace internal {

template <>
struct storage_kind_to_evaluator_kind<H2> {
  typedef H2 Kind;
};

template <>
struct storage_kind_to_shape<H2> {
typedef H2 Shape;
};

}  // end namespace internal

}  // end namespace Eigen

#endif  // __BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2ASSIGN_H__
