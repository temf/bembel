// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_H2MATRIX_H2CWISEBINARYOP_H__
#define __BEMBEL_H2MATRIX_H2CWISEBINARYOP_H__

// The contents of this file are a modification of the file
// SparseCwiseBinaryOp.h from the Eigen library
namespace Eigen {

// Here we have to handle 3 cases:
//  1 - H2 op dense
//  2 - dense op H2
//  3 - H2 op H2
// We also need to implement a 4th iterator for:
//  4 - dense op dense
// Finally, we also need to distinguish between the product and other operations
// :
//                configuration      returned mode
//  1 - H2 op dense        product      H2
//                         generic      dense
//  2 - dense op H2        product      H2
//                         generic      dense
//  3 - H2 op H2           product      H2
//                         generic      H2
//  4 - dense op dense     product      dense
//                         generic      dense

/// todo possibly adapt SparseMatrixBase to future H2MatrixBase class, once this
/// is available
template <typename BinaryOp, typename Lhs, typename Rhs>
class CwiseBinaryOpImpl<BinaryOp, Lhs, Rhs, H2>
    : public H2MatrixBase<CwiseBinaryOp<BinaryOp, Lhs, Rhs> > {
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

#if 0
// Generic "H2 OP H2"
template <typename XprType>
struct binary_H2_evaluator;

template <typename BinaryOp, typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<BinaryOp, Lhs, Rhs>, IteratorBased,
                        IteratorBased>
    : evaluator_base<CwiseBinaryOp<BinaryOp, Lhs, Rhs> > {
 protected:
  typedef typename evaluator<Lhs>::InnerIterator LhsIterator;
  typedef typename evaluator<Rhs>::InnerIterator RhsIterator;
  typedef CwiseBinaryOp<BinaryOp, Lhs, Rhs> XprType;
  typedef typename traits<XprType>::Scalar Scalar;
  typedef typename XprType::StorageIndex StorageIndex;

 public:
  class InnerIterator {
   public:
    EIGEN_STRONG_INLINE InnerIterator(const binary_evaluator& aEval,
                                      Index outer)
        : m_lhsIter(aEval.m_lhsImpl, outer),
          m_rhsIter(aEval.m_rhsImpl, outer),
          m_functor(aEval.m_functor) {
      this->operator++();
    }

    EIGEN_STRONG_INLINE InnerIterator& operator++() {
      if (m_lhsIter && m_rhsIter && (m_lhsIter.index() == m_rhsIter.index())) {
        m_id = m_lhsIter.index();
        m_value = m_functor(m_lhsIter.value(), m_rhsIter.value());
        ++m_lhsIter;
        ++m_rhsIter;
      } else if (m_lhsIter &&
                 (!m_rhsIter || (m_lhsIter.index() < m_rhsIter.index()))) {
        m_id = m_lhsIter.index();
        m_value = m_functor(m_lhsIter.value(), Scalar(0));
        ++m_lhsIter;
      } else if (m_rhsIter &&
                 (!m_lhsIter || (m_lhsIter.index() > m_rhsIter.index()))) {
        m_id = m_rhsIter.index();
        m_value = m_functor(Scalar(0), m_rhsIter.value());
        ++m_rhsIter;
      } else {
        m_value = 0;  // this is to avoid a compilation warning
        m_id = -1;
      }
      return *this;
    }

    EIGEN_STRONG_INLINE Scalar value() const { return m_value; }

    EIGEN_STRONG_INLINE StorageIndex index() const { return m_id; }
    EIGEN_STRONG_INLINE Index outer() const { return m_lhsIter.outer(); }
    EIGEN_STRONG_INLINE Index row() const {
      return Lhs::IsRowMajor ? m_lhsIter.row() : index();
    }
    EIGEN_STRONG_INLINE Index col() const {
      return Lhs::IsRowMajor ? index() : m_lhsIter.col();
    }

    EIGEN_STRONG_INLINE operator bool() const { return m_id >= 0; }

   protected:
    LhsIterator m_lhsIter;
    RhsIterator m_rhsIter;
    const BinaryOp& m_functor;
    Scalar m_value;
    StorageIndex m_id;
  };

  enum {
    CoeffReadCost = evaluator<Lhs>::CoeffReadCost +
                    evaluator<Rhs>::CoeffReadCost +
                    functor_traits<BinaryOp>::Cost,
    Flags = XprType::Flags
  };

  explicit binary_evaluator(const XprType& xpr)
      : m_functor(xpr.functor()), m_lhsImpl(xpr.lhs()), m_rhsImpl(xpr.rhs()) {
    EIGEN_INTERNAL_CHECK_COST_VALUE(functor_traits<BinaryOp>::Cost);
    EIGEN_INTERNAL_CHECK_COST_VALUE(CoeffReadCost);
  }

  inline Index nonZerosEstimate() const {
    return m_lhsImpl.nonZerosEstimate() + m_rhsImpl.nonZerosEstimate();
  }

 protected:
  const BinaryOp m_functor;
  evaluator<Lhs> m_lhsImpl;
  evaluator<Rhs> m_rhsImpl;
};

// dense op H2
template <typename BinaryOp, typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<BinaryOp, Lhs, Rhs>, IndexBased,
                        IteratorBased>
    : evaluator_base<CwiseBinaryOp<BinaryOp, Lhs, Rhs> > {
 protected:
  typedef typename evaluator<Rhs>::InnerIterator RhsIterator;
  typedef CwiseBinaryOp<BinaryOp, Lhs, Rhs> XprType;
  typedef typename traits<XprType>::Scalar Scalar;
  typedef typename XprType::StorageIndex StorageIndex;

 public:
  class InnerIterator {
    enum { IsRowMajor = (int(Rhs::Flags) & RowMajorBit) == RowMajorBit };

   public:
    EIGEN_STRONG_INLINE InnerIterator(const binary_evaluator& aEval,
                                      Index outer)
        : m_lhsEval(aEval.m_lhsImpl),
          m_rhsIter(aEval.m_rhsImpl, outer),
          m_functor(aEval.m_functor),
          m_value(0),
          m_id(-1),
          m_innerSize(aEval.m_expr.rhs().innerSize()) {
      this->operator++();
    }

    EIGEN_STRONG_INLINE InnerIterator& operator++() {
      ++m_id;
      if (m_id < m_innerSize) {
        Scalar lhsVal = m_lhsEval.coeff(IsRowMajor ? m_rhsIter.outer() : m_id,
                                        IsRowMajor ? m_id : m_rhsIter.outer());
        if (m_rhsIter && m_rhsIter.index() == m_id) {
          m_value = m_functor(lhsVal, m_rhsIter.value());
          ++m_rhsIter;
        } else
          m_value = m_functor(lhsVal, Scalar(0));
      }

      return *this;
    }

    EIGEN_STRONG_INLINE Scalar value() const {
      eigen_internal_assert(m_id < m_innerSize);
      return m_value;
    }

    EIGEN_STRONG_INLINE StorageIndex index() const { return m_id; }
    EIGEN_STRONG_INLINE Index outer() const { return m_rhsIter.outer(); }
    EIGEN_STRONG_INLINE Index row() const {
      return IsRowMajor ? m_rhsIter.outer() : m_id;
    }
    EIGEN_STRONG_INLINE Index col() const {
      return IsRowMajor ? m_id : m_rhsIter.outer();
    }

    EIGEN_STRONG_INLINE operator bool() const { return m_id < m_innerSize; }

   protected:
    const evaluator<Lhs>& m_lhsEval;
    RhsIterator m_rhsIter;
    const BinaryOp& m_functor;
    Scalar m_value;
    StorageIndex m_id;
    StorageIndex m_innerSize;
  };

  enum {
    CoeffReadCost = evaluator<Lhs>::CoeffReadCost +
                    evaluator<Rhs>::CoeffReadCost +
                    functor_traits<BinaryOp>::Cost,
    // Expose storage order of the H2 expression
    Flags = (XprType::Flags & ~RowMajorBit) | (int(Rhs::Flags) & RowMajorBit)
  };

  explicit binary_evaluator(const XprType& xpr)
      : m_functor(xpr.functor()),
        m_lhsImpl(xpr.lhs()),
        m_rhsImpl(xpr.rhs()),
        m_expr(xpr) {
    EIGEN_INTERNAL_CHECK_COST_VALUE(functor_traits<BinaryOp>::Cost);
    EIGEN_INTERNAL_CHECK_COST_VALUE(CoeffReadCost);
  }

  inline Index nonZerosEstimate() const { return m_expr.size(); }

 protected:
  const BinaryOp m_functor;
  evaluator<Lhs> m_lhsImpl;
  evaluator<Rhs> m_rhsImpl;
  const XprType& m_expr;
};

// H2 op dense
template <typename BinaryOp, typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<BinaryOp, Lhs, Rhs>, IteratorBased,
                        IndexBased>
    : evaluator_base<CwiseBinaryOp<BinaryOp, Lhs, Rhs> > {
 protected:
  typedef typename evaluator<Lhs>::InnerIterator LhsIterator;
  typedef CwiseBinaryOp<BinaryOp, Lhs, Rhs> XprType;
  typedef typename traits<XprType>::Scalar Scalar;
  typedef typename XprType::StorageIndex StorageIndex;

 public:
  class InnerIterator {
    enum { IsRowMajor = (int(Lhs::Flags) & RowMajorBit) == RowMajorBit };

   public:
    EIGEN_STRONG_INLINE InnerIterator(const binary_evaluator& aEval,
                                      Index outer)
        : m_lhsIter(aEval.m_lhsImpl, outer),
          m_rhsEval(aEval.m_rhsImpl),
          m_functor(aEval.m_functor),
          m_value(0),
          m_id(-1),
          m_innerSize(aEval.m_expr.lhs().innerSize()) {
      this->operator++();
    }

    EIGEN_STRONG_INLINE InnerIterator& operator++() {
      ++m_id;
      if (m_id < m_innerSize) {
        Scalar rhsVal = m_rhsEval.coeff(IsRowMajor ? m_lhsIter.outer() : m_id,
                                        IsRowMajor ? m_id : m_lhsIter.outer());
        if (m_lhsIter && m_lhsIter.index() == m_id) {
          m_value = m_functor(m_lhsIter.value(), rhsVal);
          ++m_lhsIter;
        } else
          m_value = m_functor(Scalar(0), rhsVal);
      }

      return *this;
    }

    EIGEN_STRONG_INLINE Scalar value() const {
      eigen_internal_assert(m_id < m_innerSize);
      return m_value;
    }

    EIGEN_STRONG_INLINE StorageIndex index() const { return m_id; }
    EIGEN_STRONG_INLINE Index outer() const { return m_lhsIter.outer(); }
    EIGEN_STRONG_INLINE Index row() const {
      return IsRowMajor ? m_lhsIter.outer() : m_id;
    }
    EIGEN_STRONG_INLINE Index col() const {
      return IsRowMajor ? m_id : m_lhsIter.outer();
    }

    EIGEN_STRONG_INLINE operator bool() const { return m_id < m_innerSize; }

   protected:
    LhsIterator m_lhsIter;
    const evaluator<Rhs>& m_rhsEval;
    const BinaryOp& m_functor;
    Scalar m_value;
    StorageIndex m_id;
    StorageIndex m_innerSize;
  };

  enum {
    CoeffReadCost = evaluator<Lhs>::CoeffReadCost +
                    evaluator<Rhs>::CoeffReadCost +
                    functor_traits<BinaryOp>::Cost,
    // Expose storage order of the H2 expression
    Flags = (XprType::Flags & ~RowMajorBit) | (int(Lhs::Flags) & RowMajorBit)
  };

  explicit binary_evaluator(const XprType& xpr)
      : m_functor(xpr.functor()),
        m_lhsImpl(xpr.lhs()),
        m_rhsImpl(xpr.rhs()),
        m_expr(xpr) {
    EIGEN_INTERNAL_CHECK_COST_VALUE(functor_traits<BinaryOp>::Cost);
    EIGEN_INTERNAL_CHECK_COST_VALUE(CoeffReadCost);
  }

  inline Index nonZerosEstimate() const { return m_expr.size(); }

 protected:
  const BinaryOp m_functor;
  evaluator<Lhs> m_lhsImpl;
  evaluator<Rhs> m_rhsImpl;
  const XprType& m_expr;
};

template <typename T,
          typename LhsKind = typename evaluator_traits<typename T::Lhs>::Kind,
          typename RhsKind = typename evaluator_traits<typename T::Rhs>::Kind,
          typename LhsScalar = typename traits<typename T::Lhs>::Scalar,
          typename RhsScalar = typename traits<typename T::Rhs>::Scalar>
struct H2_conjunction_evaluator;

// "H2 .* H2"
template <typename T1, typename T2, typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<scalar_product_op<T1, T2>, Lhs, Rhs>,
                        IteratorBased, IteratorBased>
    : H2_conjunction_evaluator<
          CwiseBinaryOp<scalar_product_op<T1, T2>, Lhs, Rhs> > {
  typedef CwiseBinaryOp<scalar_product_op<T1, T2>, Lhs, Rhs> XprType;
  typedef H2_conjunction_evaluator<XprType> Base;
  explicit binary_evaluator(const XprType& xpr) : Base(xpr) {}
};
// "dense .* H2"
template <typename T1, typename T2, typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<scalar_product_op<T1, T2>, Lhs, Rhs>,
                        IndexBased, IteratorBased>
    : H2_conjunction_evaluator<
          CwiseBinaryOp<scalar_product_op<T1, T2>, Lhs, Rhs> > {
  typedef CwiseBinaryOp<scalar_product_op<T1, T2>, Lhs, Rhs> XprType;
  typedef H2_conjunction_evaluator<XprType> Base;
  explicit binary_evaluator(const XprType& xpr) : Base(xpr) {}
};
// "H2 .* dense"
template <typename T1, typename T2, typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<scalar_product_op<T1, T2>, Lhs, Rhs>,
                        IteratorBased, IndexBased>
    : H2_conjunction_evaluator<
          CwiseBinaryOp<scalar_product_op<T1, T2>, Lhs, Rhs> > {
  typedef CwiseBinaryOp<scalar_product_op<T1, T2>, Lhs, Rhs> XprType;
  typedef H2_conjunction_evaluator<XprType> Base;
  explicit binary_evaluator(const XprType& xpr) : Base(xpr) {}
};

// "H2 ./ dense"
template <typename T1, typename T2, typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<scalar_quotient_op<T1, T2>, Lhs, Rhs>,
                        IteratorBased, IndexBased>
    : H2_conjunction_evaluator<
          CwiseBinaryOp<scalar_quotient_op<T1, T2>, Lhs, Rhs> > {
  typedef CwiseBinaryOp<scalar_quotient_op<T1, T2>, Lhs, Rhs> XprType;
  typedef H2_conjunction_evaluator<XprType> Base;
  explicit binary_evaluator(const XprType& xpr) : Base(xpr) {}
};

// "H2 && H2"
template <typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<scalar_boolean_and_op, Lhs, Rhs>,
                        IteratorBased, IteratorBased>
    : H2_conjunction_evaluator<
          CwiseBinaryOp<scalar_boolean_and_op, Lhs, Rhs> > {
  typedef CwiseBinaryOp<scalar_boolean_and_op, Lhs, Rhs> XprType;
  typedef H2_conjunction_evaluator<XprType> Base;
  explicit binary_evaluator(const XprType& xpr) : Base(xpr) {}
};
// "dense && H2"
template <typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<scalar_boolean_and_op, Lhs, Rhs>,
                        IndexBased, IteratorBased>
    : H2_conjunction_evaluator<
          CwiseBinaryOp<scalar_boolean_and_op, Lhs, Rhs> > {
  typedef CwiseBinaryOp<scalar_boolean_and_op, Lhs, Rhs> XprType;
  typedef H2_conjunction_evaluator<XprType> Base;
  explicit binary_evaluator(const XprType& xpr) : Base(xpr) {}
};
// "H2 && dense"
template <typename Lhs, typename Rhs>
struct binary_evaluator<CwiseBinaryOp<scalar_boolean_and_op, Lhs, Rhs>,
                        IteratorBased, IndexBased>
    : H2_conjunction_evaluator<
          CwiseBinaryOp<scalar_boolean_and_op, Lhs, Rhs> > {
  typedef CwiseBinaryOp<scalar_boolean_and_op, Lhs, Rhs> XprType;
  typedef H2_conjunction_evaluator<XprType> Base;
  explicit binary_evaluator(const XprType& xpr) : Base(xpr) {}
};

// "H2 ^ H2"
template <typename XprType>
struct H2_conjunction_evaluator<XprType, IteratorBased, IteratorBased>
    : evaluator_base<XprType> {
 protected:
  typedef typename XprType::Functor BinaryOp;
  typedef typename XprType::Lhs LhsArg;
  typedef typename XprType::Rhs RhsArg;
  typedef typename evaluator<LhsArg>::InnerIterator LhsIterator;
  typedef typename evaluator<RhsArg>::InnerIterator RhsIterator;
  typedef typename XprType::StorageIndex StorageIndex;
  typedef typename traits<XprType>::Scalar Scalar;

 public:
  class InnerIterator {
   public:
    EIGEN_STRONG_INLINE InnerIterator(const H2_conjunction_evaluator& aEval,
                                      Index outer)
        : m_lhsIter(aEval.m_lhsImpl, outer),
          m_rhsIter(aEval.m_rhsImpl, outer),
          m_functor(aEval.m_functor) {
      while (m_lhsIter && m_rhsIter &&
             (m_lhsIter.index() != m_rhsIter.index())) {
        if (m_lhsIter.index() < m_rhsIter.index())
          ++m_lhsIter;
        else
          ++m_rhsIter;
      }
    }

    EIGEN_STRONG_INLINE InnerIterator& operator++() {
      ++m_lhsIter;
      ++m_rhsIter;
      while (m_lhsIter && m_rhsIter &&
             (m_lhsIter.index() != m_rhsIter.index())) {
        if (m_lhsIter.index() < m_rhsIter.index())
          ++m_lhsIter;
        else
          ++m_rhsIter;
      }
      return *this;
    }

    EIGEN_STRONG_INLINE Scalar value() const {
      return m_functor(m_lhsIter.value(), m_rhsIter.value());
    }

    EIGEN_STRONG_INLINE StorageIndex index() const { return m_lhsIter.index(); }
    EIGEN_STRONG_INLINE Index outer() const { return m_lhsIter.outer(); }
    EIGEN_STRONG_INLINE Index row() const { return m_lhsIter.row(); }
    EIGEN_STRONG_INLINE Index col() const { return m_lhsIter.col(); }

    EIGEN_STRONG_INLINE operator bool() const {
      return (m_lhsIter && m_rhsIter);
    }

   protected:
    LhsIterator m_lhsIter;
    RhsIterator m_rhsIter;
    const BinaryOp& m_functor;
  };

  enum {
    CoeffReadCost = evaluator<LhsArg>::CoeffReadCost +
                    evaluator<RhsArg>::CoeffReadCost +
                    functor_traits<BinaryOp>::Cost,
    Flags = XprType::Flags
  };

  explicit H2_conjunction_evaluator(const XprType& xpr)
      : m_functor(xpr.functor()), m_lhsImpl(xpr.lhs()), m_rhsImpl(xpr.rhs()) {
    EIGEN_INTERNAL_CHECK_COST_VALUE(functor_traits<BinaryOp>::Cost);
    EIGEN_INTERNAL_CHECK_COST_VALUE(CoeffReadCost);
  }

  inline Index nonZerosEstimate() const {
    return (std::min)(m_lhsImpl.nonZerosEstimate(),
                      m_rhsImpl.nonZerosEstimate());
  }

 protected:
  const BinaryOp m_functor;
  evaluator<LhsArg> m_lhsImpl;
  evaluator<RhsArg> m_rhsImpl;
};

// "dense ^ H2"
template <typename XprType>
struct H2_conjunction_evaluator<XprType, IndexBased, IteratorBased>
    : evaluator_base<XprType> {
 protected:
  typedef typename XprType::Functor BinaryOp;
  typedef typename XprType::Lhs LhsArg;
  typedef typename XprType::Rhs RhsArg;
  typedef evaluator<LhsArg> LhsEvaluator;
  typedef typename evaluator<RhsArg>::InnerIterator RhsIterator;
  typedef typename XprType::StorageIndex StorageIndex;
  typedef typename traits<XprType>::Scalar Scalar;

 public:
  class InnerIterator {
    enum { IsRowMajor = (int(RhsArg::Flags) & RowMajorBit) == RowMajorBit };

   public:
    EIGEN_STRONG_INLINE InnerIterator(const H2_conjunction_evaluator& aEval,
                                      Index outer)
        : m_lhsEval(aEval.m_lhsImpl),
          m_rhsIter(aEval.m_rhsImpl, outer),
          m_functor(aEval.m_functor),
          m_outer(outer) {}

    EIGEN_STRONG_INLINE InnerIterator& operator++() {
      ++m_rhsIter;
      return *this;
    }

    EIGEN_STRONG_INLINE Scalar value() const {
      return m_functor(
          m_lhsEval.coeff(IsRowMajor ? m_outer : m_rhsIter.index(),
                          IsRowMajor ? m_rhsIter.index() : m_outer),
          m_rhsIter.value());
    }

    EIGEN_STRONG_INLINE StorageIndex index() const { return m_rhsIter.index(); }
    EIGEN_STRONG_INLINE Index outer() const { return m_rhsIter.outer(); }
    EIGEN_STRONG_INLINE Index row() const { return m_rhsIter.row(); }
    EIGEN_STRONG_INLINE Index col() const { return m_rhsIter.col(); }

    EIGEN_STRONG_INLINE operator bool() const { return m_rhsIter; }

   protected:
    const LhsEvaluator& m_lhsEval;
    RhsIterator m_rhsIter;
    const BinaryOp& m_functor;
    const Index m_outer;
  };

  enum {
    CoeffReadCost = evaluator<LhsArg>::CoeffReadCost +
                    evaluator<RhsArg>::CoeffReadCost +
                    functor_traits<BinaryOp>::Cost,
    // Expose storage order of the H2 expression
    Flags = (XprType::Flags & ~RowMajorBit) | (int(RhsArg::Flags) & RowMajorBit)
  };

  explicit H2_conjunction_evaluator(const XprType& xpr)
      : m_functor(xpr.functor()), m_lhsImpl(xpr.lhs()), m_rhsImpl(xpr.rhs()) {
    EIGEN_INTERNAL_CHECK_COST_VALUE(functor_traits<BinaryOp>::Cost);
    EIGEN_INTERNAL_CHECK_COST_VALUE(CoeffReadCost);
  }

  inline Index nonZerosEstimate() const { return m_rhsImpl.nonZerosEstimate(); }

 protected:
  const BinaryOp m_functor;
  evaluator<LhsArg> m_lhsImpl;
  evaluator<RhsArg> m_rhsImpl;
};

// "H2 ^ dense"
template <typename XprType>
struct H2_conjunction_evaluator<XprType, IteratorBased, IndexBased>
    : evaluator_base<XprType> {
 protected:
  typedef typename XprType::Functor BinaryOp;
  typedef typename XprType::Lhs LhsArg;
  typedef typename XprType::Rhs RhsArg;
  typedef typename evaluator<LhsArg>::InnerIterator LhsIterator;
  typedef evaluator<RhsArg> RhsEvaluator;
  typedef typename XprType::StorageIndex StorageIndex;
  typedef typename traits<XprType>::Scalar Scalar;

 public:
  class InnerIterator {
    enum { IsRowMajor = (int(LhsArg::Flags) & RowMajorBit) == RowMajorBit };

   public:
    EIGEN_STRONG_INLINE InnerIterator(const H2_conjunction_evaluator& aEval,
                                      Index outer)
        : m_lhsIter(aEval.m_lhsImpl, outer),
          m_rhsEval(aEval.m_rhsImpl),
          m_functor(aEval.m_functor),
          m_outer(outer) {}

    EIGEN_STRONG_INLINE InnerIterator& operator++() {
      ++m_lhsIter;
      return *this;
    }

    EIGEN_STRONG_INLINE Scalar value() const {
      return m_functor(
          m_lhsIter.value(),
          m_rhsEval.coeff(IsRowMajor ? m_outer : m_lhsIter.index(),
                          IsRowMajor ? m_lhsIter.index() : m_outer));
    }

    EIGEN_STRONG_INLINE StorageIndex index() const { return m_lhsIter.index(); }
    EIGEN_STRONG_INLINE Index outer() const { return m_lhsIter.outer(); }
    EIGEN_STRONG_INLINE Index row() const { return m_lhsIter.row(); }
    EIGEN_STRONG_INLINE Index col() const { return m_lhsIter.col(); }

    EIGEN_STRONG_INLINE operator bool() const { return m_lhsIter; }

   protected:
    LhsIterator m_lhsIter;
    const evaluator<RhsArg>& m_rhsEval;
    const BinaryOp& m_functor;
    const Index m_outer;
  };

  enum {
    CoeffReadCost = evaluator<LhsArg>::CoeffReadCost +
                    evaluator<RhsArg>::CoeffReadCost +
                    functor_traits<BinaryOp>::Cost,
    // Expose storage order of the H2 expression
    Flags = (XprType::Flags & ~RowMajorBit) | (int(LhsArg::Flags) & RowMajorBit)
  };

  explicit H2_conjunction_evaluator(const XprType& xpr)
      : m_functor(xpr.functor()), m_lhsImpl(xpr.lhs()), m_rhsImpl(xpr.rhs()) {
    EIGEN_INTERNAL_CHECK_COST_VALUE(functor_traits<BinaryOp>::Cost);
    EIGEN_INTERNAL_CHECK_COST_VALUE(CoeffReadCost);
  }

  inline Index nonZerosEstimate() const { return m_lhsImpl.nonZerosEstimate(); }

 protected:
  const BinaryOp m_functor;
  evaluator<LhsArg> m_lhsImpl;
  evaluator<RhsArg> m_rhsImpl;
};
#endif

}  // namespace internal

/***************************************************************************
 * Implementation of H2MatrixBase and H2Cwise functions/operators
 ***************************************************************************/

#if 0
template <typename Derived>
template <typename OtherDerived>
Derived& H2MatrixBase<Derived>::operator+=(
    const EigenBase<OtherDerived>& other) {
  call_assignment(
      derived(), other.derived(),
      internal::add_assign_op<Scalar, typename OtherDerived::Scalar>());
  return derived();
}

template <typename Derived>
template <typename OtherDerived>
Derived& H2MatrixBase<Derived>::operator-=(
    const EigenBase<OtherDerived>& other) {
  call_assignment(derived(), other.derived(),
                  internal::assign_op<Scalar, typename OtherDerived::Scalar>());
  return derived();
}

template <typename Derived>
template <typename OtherDerived>
EIGEN_STRONG_INLINE Derived& H2MatrixBase<Derived>::operator-=(
    const H2MatrixBase<OtherDerived>& other) {
  return derived() = derived() - other.derived();
}

template <typename Derived>
template <typename OtherDerived>
EIGEN_STRONG_INLINE Derived& H2MatrixBase<Derived>::operator+=(
    const H2MatrixBase<OtherDerived>& other) {
  return derived() = derived() + other.derived();
}

template <typename Derived>
template <typename OtherDerived>
Derived& H2MatrixBase<Derived>::operator+=(
    const DiagonalBase<OtherDerived>& other) {
  call_assignment_no_alias(
      derived(), other.derived(),
      internal::add_assign_op<Scalar, typename OtherDerived::Scalar>());
  return derived();
}

template <typename Derived>
template <typename OtherDerived>
Derived& H2MatrixBase<Derived>::operator-=(
    const DiagonalBase<OtherDerived>& other) {
  call_assignment_no_alias(
      derived(), other.derived(),
      internal::sub_assign_op<Scalar, typename OtherDerived::Scalar>());
  return derived();
}

template <typename Derived>
template <typename OtherDerived>
EIGEN_STRONG_INLINE const typename H2MatrixBase<
    Derived>::template CwiseProductDenseReturnType<OtherDerived>::Type
H2MatrixBase<Derived>::cwiseProduct(
    const MatrixBase<OtherDerived>& other) const {
  return typename CwiseProductDenseReturnType<OtherDerived>::Type(
      derived(), other.derived());
}
#endif

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

}  // end namespace Eigen

#endif
