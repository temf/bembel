// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2UTIL_H__
#define __BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2UTIL_H__

// The contents of this file are a modification of the file
// SparseUtil.h from the Eigen library
namespace Eigen {

#if 0

#ifdef NDEBUG
#define EIGEN_DBG_H2(X)
#else
#define EIGEN_DBG_H2(X) X
#endif

#define EIGEN_H2_INHERIT_ASSIGNMENT_OPERATOR(Derived, Op)          \
  template <typename OtherDerived>                                 \
  EIGEN_STRONG_INLINE Derived& operator Op(                        \
      const Eigen::H2MatrixBase<OtherDerived>& other) {            \
    return Base::operator Op(other.derived());                     \
  }                                                                \
  EIGEN_STRONG_INLINE Derived& operator Op(const Derived& other) { \
    return Base::operator Op(other);                               \
  }

#define EIGEN_H2_INHERIT_SCALAR_ASSIGNMENT_OPERATOR(Derived, Op)  \
  template <typename Other>                                       \
  EIGEN_STRONG_INLINE Derived& operator Op(const Other& scalar) { \
    return Base::operator Op(scalar);                             \
  }

#define EIGEN_H2_INHERIT_ASSIGNMENT_OPERATORS(Derived) \
  EIGEN_H2_INHERIT_ASSIGNMENT_OPERATOR(Derived, =)

#define EIGEN_H2_PUBLIC_INTERFACE(Derived) \
  EIGEN_GENERIC_PUBLIC_INTERFACE(Derived)

const int CoherentAccessPattern = 0x1;
const int InnerRandomAccessPattern = 0x2 | CoherentAccessPattern;
const int OuterRandomAccessPattern = 0x4 | CoherentAccessPattern;
const int RandomAccessPattern =
    0x8 | OuterRandomAccessPattern | InnerRandomAccessPattern;

template <typename _Scalar, int _Flags = 0, typename _StorageIndex = int>
class H2Matrix;
template <typename _Scalar, int _Flags = 0, typename _StorageIndex = int>
class DynamicH2Matrix;
template <typename _Scalar, int _Flags = 0, typename _StorageIndex = int>
class H2Vector;
template <typename _Scalar, int _Flags = 0, typename _StorageIndex = int>
class MappedH2Matrix;

template <typename MatrixType, unsigned int UpLo>
class H2SelfAdjointView;
template <typename Lhs, typename Rhs>
class H2DiagonalProduct;
template <typename MatrixType>
class H2View;

template <typename Lhs, typename Rhs>
class H2H2Product;
template <typename Lhs, typename Rhs>
class H2TimeDenseProduct;
template <typename Lhs, typename Rhs>
class DenseTimeH2Product;
template <typename Lhs, typename Rhs, bool Transpose>
class H2DenseOuterProduct;

template <typename Lhs, typename Rhs>
struct H2H2ProductReturnType;
template <typename Lhs, typename Rhs,
          int InnerSize = EIGEN_SIZE_MIN_PREFER_FIXED(
              internal::traits<Lhs>::ColsAtCompileTime,
              internal::traits<Rhs>::RowsAtCompileTime)>
struct DenseH2ProductReturnType;

template <typename Lhs, typename Rhs,
          int InnerSize = EIGEN_SIZE_MIN_PREFER_FIXED(
              internal::traits<Lhs>::ColsAtCompileTime,
              internal::traits<Rhs>::RowsAtCompileTime)>
struct H2DenseProductReturnType;
template <typename MatrixType, int UpLo>
class H2SymmetricPermutationProduct;

#endif

namespace internal {

template <typename T, int Rows, int Cols, int Flags>
struct H2_eval;

#if 0
template <typename T>
struct eval<T, H2>
    : H2_eval<T, traits<T>::RowsAtCompileTime, traits<T>::ColsAtCompileTime,
                  traits<T>::Flags> {};

template <typename T, int Cols, int Flags>
struct H2_eval<T, 1, Cols, Flags> {
  typedef typename traits<T>::Scalar _Scalar;
  typedef typename traits<T>::StorageIndex _StorageIndex;

 public:
  typedef H2Vector<_Scalar, RowMajor, _StorageIndex> type;
};

template <typename T, int Rows, int Flags>
struct H2_eval<T, Rows, 1, Flags> {
  typedef typename traits<T>::Scalar _Scalar;
  typedef typename traits<T>::StorageIndex _StorageIndex;

 public:
  typedef H2Vector<_Scalar, ColMajor, _StorageIndex> type;
};
#endif

// TODO this seems almost identical to plain_matrix_type<T, H2>
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

#if 0
template <typename T, int Flags>
struct H2_eval<T, 1, 1, Flags> {
  typedef typename traits<T>::Scalar _Scalar;

 public:
  typedef Matrix<_Scalar, 1, 1> type;
};

template <typename T>
struct plain_matrix_type<T, H2> {
  typedef typename traits<T>::Scalar _Scalar;
  typedef typename traits<T>::StorageIndex _StorageIndex;
  enum {
    _Options = ((evaluator<T>::Flags & RowMajorBit) == RowMajorBit) ? RowMajor
                                                                    : ColMajor
  };

 public:
  typedef H2Matrix<_Scalar, _Options, _StorageIndex> type;
};
#endif

template <typename T>
struct plain_object_eval<T, H2>
    : H2_eval<T, traits<T>::RowsAtCompileTime, traits<T>::ColsAtCompileTime,
              evaluator<T>::Flags> {};

#if 0
template <typename Decomposition, typename RhsType>
struct solve_traits<Decomposition, RhsType, H2> {
  typedef typename H2_eval<RhsType, RhsType::RowsAtCompileTime,
                               RhsType::ColsAtCompileTime,
                               traits<RhsType>::Flags>::type PlainObject;
};

template <typename Derived>
struct generic_xpr_base<Derived, MatrixXpr, H2> {
  typedef H2MatrixBase<Derived> type;
};

struct H2TriangularShape {
  static std::string debugName() { return "H2TriangularShape"; }
};
struct H2SelfAdjointShape {
  static std::string debugName() { return "H2SelfAdjointShape"; }
};

template <>
struct glue_shapes<H2Shape, SelfAdjointShape> {
  typedef H2SelfAdjointShape type;
};
template <>
struct glue_shapes<H2Shape, TriangularShape> {
  typedef H2TriangularShape type;
};

#endif

}  // end namespace internal

}  // end namespace Eigen

#endif  // EIGEN_H2UTIL_H
