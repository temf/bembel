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

#if 0
template <typename Derived>
template <typename OtherDerived>
Derived &H2MatrixBase<Derived>::operator=(
    const EigenBase<OtherDerived> &other) {
  internal::call_assignment_no_alias(derived(), other.derived());
  return derived();
}

template <typename Derived>
template <typename OtherDerived>
Derived &H2MatrixBase<Derived>::operator=(
    const ReturnByValue<OtherDerived> &other) {
  // TODO use the evaluator mechanism
  other.evalTo(derived());
  return derived();
}

template <typename Derived>
template <typename OtherDerived>
inline Derived &H2MatrixBase<Derived>::operator=(
    const H2MatrixBase<OtherDerived> &other) {
  // by default sparse evaluation do not alias, so we can safely bypass the
  // generic call_assignment routine
  internal::Assignment<
      Derived, OtherDerived,
      internal::assign_op<Scalar, typename OtherDerived::Scalar> >::
      run(derived(), other.derived(),
          internal::assign_op<Scalar, typename OtherDerived::Scalar>());
  return derived();
}

template <typename Derived>
inline Derived &H2MatrixBase<Derived>::operator=(const Derived &other) {
  internal::call_assignment_no_alias(derived(), other.derived());
  return derived();
}
#endif

namespace internal {

#if 0
template <>
struct storage_kind_to_evaluator_kind<H2> {
  typedef IteratorBased Kind;
};
#endif

template <>
struct storage_kind_to_shape<H2> {
typedef H2 Shape;
};

#if 0
struct H22H2 {};
struct H22Dense {};

template <>
struct AssignmentKind<H2, H2> {
  typedef H22H2 Kind;
};
template <>
struct AssignmentKind<H2, H2TriangularShape> {
  typedef H22H2 Kind;
};
template <>
struct AssignmentKind<DenseShape, H2> {
  typedef H22Dense Kind;
};
template <>
struct AssignmentKind<DenseShape, H2TriangularShape> {
  typedef H22Dense Kind;
};

template <typename DstXprType, typename SrcXprType>
void assign_H2_to_H2(DstXprType &dst, const SrcXprType &src) {
  // implementation to be done
}

// Generic H2 to H2 assignment
template <typename DstXprType, typename SrcXprType, typename Functor>
struct Assignment<DstXprType, SrcXprType, Functor, H22H2> {
  static void run(
      DstXprType &dst, const SrcXprType &src,
      const internal::assign_op<typename DstXprType::Scalar,
                                typename SrcXprType::Scalar> & /*func*/) {
    assign_H2_to_H2(dst.derived(), src.derived());
  }
};

// Generic H2 to Dense assignment
template <typename DstXprType, typename SrcXprType, typename Functor>
struct Assignment<DstXprType, SrcXprType, Functor, H22Dense> {
  static void run(DstXprType &dst, const SrcXprType &src, const Functor &func) {
    // implementation to be done
  }
};

// Specialization for "dst = dec.solve(rhs)"
// NOTE we need to specialize it for H22H2 to avoid ambiguous
// specialization error
template <typename DstXprType, typename DecType, typename RhsType,
          typename Scalar>
struct Assignment<DstXprType, Solve<DecType, RhsType>,
                  internal::assign_op<Scalar, Scalar>, H22H2> {
  typedef Solve<DecType, RhsType> SrcXprType;
  static void run(DstXprType &dst, const SrcXprType &src,
                  const internal::assign_op<Scalar, Scalar> &) {
    Index dstRows = src.rows();
    Index dstCols = src.cols();
    if ((dst.rows() != dstRows) || (dst.cols() != dstCols))
      dst.resize(dstRows, dstCols);

    src.dec()._solve_impl(src.rhs(), dst);
  }
};

struct Diagonal2H2 {};

template <>
struct AssignmentKind<H2, DiagonalShape> {
  typedef Diagonal2H2 Kind;
};

template <typename DstXprType, typename SrcXprType, typename Functor>
struct Assignment<DstXprType, SrcXprType, Functor, Diagonal2H2> {
  typedef typename DstXprType::StorageIndex StorageIndex;
  typedef typename DstXprType::Scalar Scalar;
  typedef Array<StorageIndex, Dynamic, 1> ArrayXI;
  typedef Array<Scalar, Dynamic, 1> ArrayXS;
  template <int Options>
  static void run(
      H2Matrix<Scalar, Options, StorageIndex> &dst, const SrcXprType &src,
      const internal::assign_op<typename DstXprType::Scalar,
                                typename SrcXprType::Scalar> & /*func*/) {
    // implementation to be done
  }

  template <typename DstDerived>
  static void run(
      H2MatrixBase<DstDerived> &dst, const SrcXprType &src,
      const internal::assign_op<typename DstXprType::Scalar,
                                typename SrcXprType::Scalar> & /*func*/) {
    dst.diagonal() = src.diagonal();
  }

  static void run(
      DstXprType &dst, const SrcXprType &src,
      const internal::add_assign_op<typename DstXprType::Scalar,
                                    typename SrcXprType::Scalar> & /*func*/) {
    dst.diagonal() += src.diagonal();
  }

  static void run(
      DstXprType &dst, const SrcXprType &src,
      const internal::sub_assign_op<typename DstXprType::Scalar,
                                    typename SrcXprType::Scalar> & /*func*/) {
    dst.diagonal() -= src.diagonal();
  }
};
#endif
}  // end namespace internal

}  // end namespace Eigen

#endif  // __BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2ASSIGN_H__
