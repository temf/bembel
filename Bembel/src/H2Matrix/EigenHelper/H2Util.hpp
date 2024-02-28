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
#ifndef BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2UTIL_HPP_
#define BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2UTIL_HPP_

// The contents of this file are a modification of the file
// SparseUtil.h from the Eigen library
namespace Eigen {

namespace internal {

template <typename T, int Rows, int Cols, int Flags>
struct H2_eval;

template <typename T, int Rows, int Cols, int Flags>
struct H2_eval {
  typedef typename traits<T>::Scalar _Scalar;
  typedef typename traits<T>::StorageIndex _StorageIndex;
  enum {
    _Options = ((Flags & RowMajorBit) == RowMajorBit) ? RowMajor : ColMajor
  };

 public:
  typedef H2Matrix<_Scalar> type;
};

template <typename T>
struct plain_object_eval<T, H2>
    : H2_eval<T, traits<T>::RowsAtCompileTime, traits<T>::ColsAtCompileTime,
              evaluator<T>::Flags> {};

}  // end namespace internal

}  // end namespace Eigen

#endif  // BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2UTIL_HPP_
