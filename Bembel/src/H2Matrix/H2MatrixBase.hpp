// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_H2MATRIX_H2MATRIXBASE_H__
#define __BEMBEL_H2MATRIX_H2MATRIXBASE_H__

namespace Eigen {

template <typename Derived>
class H2Matrix;

/** \ingroup H2Matrix
 */
template <typename Derived>
class H2MatrixBase : public EigenBase<Derived> {
 public:
  typedef typename internal::traits<Derived>::Scalar Scalar;
  typedef typename internal::ref_selector<Derived>::type Nested;
  typedef Scalar CoeffReturnType;

  enum {
    RowsAtCompileTime = internal::traits<Derived>::RowsAtCompileTime,
    ColsAtCompileTime = internal::traits<Derived>::ColsAtCompileTime,
    SizeAtCompileTime = (internal::size_at_compile_time<
                         internal::traits<Derived>::RowsAtCompileTime,
                         internal::traits<Derived>::ColsAtCompileTime>::ret),
    MaxRowsAtCompileTime = RowsAtCompileTime,
    MaxColsAtCompileTime = ColsAtCompileTime,
    MaxSizeAtCompileTime =
        (internal::size_at_compile_time<MaxRowsAtCompileTime,
                                        MaxColsAtCompileTime>::ret),
    IsVectorAtCompileTime = false
  };

  inline const Derived& derived() const {
    return *static_cast<const Derived*>(this);
  }
  inline Derived& derived() { return *static_cast<Derived*>(this); }

  // H2 * dense
  template <typename OtherDerived>
  const Product<Derived, OtherDerived> operator*(
      const MatrixBase<OtherDerived>& other) const {
    return Product<Derived, OtherDerived>(derived(), other.derived());
  }

 private:
  // these operations cannot be performed exactly, so we declare them private
  template <typename OtherDerived>
  Derived& operator+=(const SparseMatrixBase<OtherDerived>& other);
  template <typename OtherDerived>
  Derived& operator-=(const SparseMatrixBase<OtherDerived>& other);
  template <typename OtherDerived>
  Derived& operator+=(const DiagonalBase<OtherDerived>& other);
  template <typename OtherDerived>
  Derived& operator-=(const DiagonalBase<OtherDerived>& other);
  template <typename OtherDerived>
  Derived& operator+=(const EigenBase<OtherDerived>& other);
  template <typename OtherDerived>
  Derived& operator-=(const EigenBase<OtherDerived>& other);

  Derived& operator*=(const Scalar& other);
  Derived& operator/=(const Scalar& other);
};

namespace internal {

// adaption from SparseMatrixBase from Eigen
template <typename Derived>
struct evaluator<H2MatrixBase<Derived>> : evaluator_base<Derived> {
  typedef typename Derived::Scalar Scalar;

  enum { CoeffReadCost = NumTraits<Scalar>::ReadCost, Flags = Derived::Flags };

  operator Derived&() { return m_matrix->const_cast_derived(); }
  operator const Derived&() const { return *m_matrix; }

  typedef typename DenseCoeffsBase<Derived, ReadOnlyAccessors>::CoeffReturnType
      CoeffReturnType;

  const Derived* m_matrix;
  const Scalar m_zero;
};

}  // end namespace internal

}  // end namespace Eigen

#endif
