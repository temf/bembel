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
#ifndef BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2CWISEBINARYOP_HPP_
#define BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2CWISEBINARYOP_HPP_

// The contents of this file are a modification of the file
// SparseCwiseBinaryOp.h from the Eigen library
namespace Eigen {

template <typename BinaryOp, typename Lhs, typename Rhs>
class CwiseBinaryOpImpl<BinaryOp, Lhs, Rhs, H2>
    : public H2MatrixBase<CwiseBinaryOp<BinaryOp, Lhs, Rhs>> {
 public:
  typedef CwiseBinaryOp<BinaryOp, Lhs, Rhs> Derived;
  typedef H2MatrixBase<Derived> Base;
  EIGEN_GENERIC_PUBLIC_INTERFACE(Derived)
  CwiseBinaryOpImpl() {
    EIGEN_STATIC_ASSERT(
        ((!internal::is_same<
             typename internal::traits<Lhs>::StorageKind,
             typename internal::traits<Rhs>::StorageKind>::value) ||
         ((internal::evaluator<Lhs>::Flags & RowMajorBit) ==
          (internal::evaluator<Rhs>::Flags & RowMajorBit))),
        THE_STORAGE_ORDER_OF_BOTH_SIDES_MUST_MATCH);
  }
};

namespace internal {

// brute force solution to allow the use of CwiseBinaryOps in solvers
template <typename BinaryOp, typename Lhs, typename Rhs>
struct is_ref_compatible<CwiseBinaryOp<BinaryOp, Lhs, Rhs>> {
  enum { value = false };
};
template <typename BinaryOp, typename Lhs, typename Rhs>
struct is_ref_compatible<const CwiseBinaryOp<BinaryOp, Lhs, Rhs>>
    : is_ref_compatible<CwiseBinaryOp<BinaryOp, Lhs, Rhs>> {};

// allow H2+H2, H2-H2, etc.
template <typename BinaryOp, typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<BinaryOp, Lhs, Rhs>, H2, H2>
    : evaluator_base<CwiseBinaryOp<BinaryOp, Lhs, Rhs>> {
 protected:
  typedef CwiseBinaryOp<BinaryOp, Lhs, Rhs> XprType;
  typedef typename traits<XprType>::Scalar Scalar;
  typedef typename XprType::StorageIndex StorageIndex;

 public:
  enum {
    CoeffReadCost = HugeCost,
    Flags = (XprType::Flags & ~RowMajorBit) | (int(Rhs::Flags) & RowMajorBit)
  };

  explicit binary_evaluator(const XprType& xpr)
      : m_functor(xpr.functor()),
        m_lhsImpl(xpr.lhs()),
        m_rhsImpl(xpr.rhs()),
        m_expr(xpr) {
    EIGEN_INTERNAL_CHECK_COST_VALUE(CoeffReadCost);
  }

 protected:
  const BinaryOp m_functor;
  evaluator<Lhs> m_lhsImpl;
  evaluator<Rhs> m_rhsImpl;
  const XprType& m_expr;
};
// interaction structs with dense
template <typename BinaryOp, typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<BinaryOp, Lhs, Rhs>, IndexBased, H2>
    : binary_evaluator<CwiseBinaryOp<BinaryOp, Lhs, Rhs>, H2, H2> {};
template <typename BinaryOp, typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<BinaryOp, Lhs, Rhs>, H2, IndexBased>
    : binary_evaluator<CwiseBinaryOp<BinaryOp, Lhs, Rhs>, H2, H2> {};
// interaction structs with sparse
template <typename BinaryOp, typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<BinaryOp, Lhs, Rhs>, IteratorBased, H2>
    : binary_evaluator<CwiseBinaryOp<BinaryOp, Lhs, Rhs>, H2, H2> {};
template <typename BinaryOp, typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<BinaryOp, Lhs, Rhs>, H2, IteratorBased>
    : binary_evaluator<CwiseBinaryOp<BinaryOp, Lhs, Rhs>, H2, H2> {};

}  // namespace internal

/***************************************************************************
 * Implementation of H2MatrixBase and H2Cwise functions/operators
 ***************************************************************************/

template <typename DenseDerived, typename H2Derived>
EIGEN_STRONG_INLINE const
    CwiseBinaryOp<internal::scalar_sum_op<typename DenseDerived::Scalar,
                                          typename H2Derived::Scalar>,
                  const DenseDerived, const H2Derived>
    operator+(const MatrixBase<DenseDerived>& a,
              const H2MatrixBase<H2Derived>& b) {
  return CwiseBinaryOp<internal::scalar_sum_op<typename DenseDerived::Scalar,
                                               typename H2Derived::Scalar>,
                       const DenseDerived, const H2Derived>(a.derived(),
                                                            b.derived());
}

template <typename H2Derived, typename DenseDerived>
EIGEN_STRONG_INLINE const
    CwiseBinaryOp<internal::scalar_sum_op<typename H2Derived::Scalar,
                                          typename DenseDerived::Scalar>,
                  const H2Derived, const DenseDerived>
    operator+(const H2MatrixBase<H2Derived>& a,
              const MatrixBase<DenseDerived>& b) {
  return CwiseBinaryOp<internal::scalar_sum_op<typename H2Derived::Scalar,
                                               typename DenseDerived::Scalar>,
                       const H2Derived, const DenseDerived>(a.derived(),
                                                            b.derived());
}

template <typename SparseDerived, typename H2Derived>
EIGEN_STRONG_INLINE const
    CwiseBinaryOp<internal::scalar_sum_op<typename SparseDerived::Scalar,
                                          typename H2Derived::Scalar>,
                  const SparseDerived, const H2Derived>
    operator+(const SparseMatrixBase<SparseDerived>& a,
              const H2MatrixBase<H2Derived>& b) {
  return CwiseBinaryOp<internal::scalar_sum_op<typename SparseDerived::Scalar,
                                               typename H2Derived::Scalar>,
                       const SparseDerived, const H2Derived>(a.derived(),
                                                             b.derived());
}

template <typename H2Derived, typename SparseDerived>
EIGEN_STRONG_INLINE const
    CwiseBinaryOp<internal::scalar_sum_op<typename H2Derived::Scalar,
                                          typename SparseDerived::Scalar>,
                  const H2Derived, const SparseDerived>
    operator+(const H2MatrixBase<H2Derived>& a,
              const SparseMatrixBase<SparseDerived>& b) {
  return CwiseBinaryOp<internal::scalar_sum_op<typename H2Derived::Scalar,
                                               typename SparseDerived::Scalar>,
                       const H2Derived, const SparseDerived>(a.derived(),
                                                             b.derived());
}

template <typename DenseDerived, typename H2Derived>
EIGEN_STRONG_INLINE const
    CwiseBinaryOp<internal::scalar_difference_op<typename DenseDerived::Scalar,
                                                 typename H2Derived::Scalar>,
                  const DenseDerived, const H2Derived>
    operator-(const MatrixBase<DenseDerived>& a,
              const H2MatrixBase<H2Derived>& b) {
  return CwiseBinaryOp<
      internal::scalar_difference_op<typename DenseDerived::Scalar,
                                     typename H2Derived::Scalar>,
      const DenseDerived, const H2Derived>(a.derived(), b.derived());
}

template <typename H2Derived, typename DenseDerived>
EIGEN_STRONG_INLINE const
    CwiseBinaryOp<internal::scalar_difference_op<typename H2Derived::Scalar,
                                                 typename DenseDerived::Scalar>,
                  const H2Derived, const DenseDerived>
    operator-(const H2MatrixBase<H2Derived>& a,
              const MatrixBase<DenseDerived>& b) {
  return CwiseBinaryOp<
      internal::scalar_difference_op<typename H2Derived::Scalar,
                                     typename DenseDerived::Scalar>,
      const H2Derived, const DenseDerived>(a.derived(), b.derived());
}

template <typename SparseDerived, typename H2Derived>
EIGEN_STRONG_INLINE const
    CwiseBinaryOp<internal::scalar_difference_op<typename SparseDerived::Scalar,
                                                 typename H2Derived::Scalar>,
                  const SparseDerived, const H2Derived>
    operator-(const SparseMatrixBase<SparseDerived>& a,
              const H2MatrixBase<H2Derived>& b) {
  return CwiseBinaryOp<
      internal::scalar_difference_op<typename SparseDerived::Scalar,
                                     typename H2Derived::Scalar>,
      const SparseDerived, const H2Derived>(a.derived(), b.derived());
}

template <typename H2Derived, typename SparseDerived>
EIGEN_STRONG_INLINE const CwiseBinaryOp<
    internal::scalar_difference_op<typename H2Derived::Scalar,
                                   typename SparseDerived::Scalar>,
    const H2Derived, const SparseDerived>
operator-(const H2MatrixBase<H2Derived>& a,
          const SparseMatrixBase<SparseDerived>& b) {
  return CwiseBinaryOp<
      internal::scalar_difference_op<typename H2Derived::Scalar,
                                     typename SparseDerived::Scalar>,
      const H2Derived, const SparseDerived>(a.derived(), b.derived());
}

}  // end namespace Eigen

#endif  // BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2CWISEBINARYOP_HPP_
