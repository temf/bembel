// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_HIERARCHICALMATRIX__
#define __BEMBEL_HIERARCHICALMATRIX__

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "Discretization.hpp"
#include "cluster_tree.h"
#include "hmatrixfactory.h"
#include "hmatrixsettings.h"

/**
 *  \namespace DirtyLittleHelpers
 *  \brief Provides useful stuff like compile time checks for types.
 *         Currently, this is not used.
 **/
namespace DirtyLittleHelpers {
template <bool condition>
struct testMatchingTypes {};
template <>
struct testMatchingTypes<true> {
  enum {
    YOU_CHOSE_INCOMPATIBLE_TYPES_FOR_HIERARCHICALMATRIX_AND_PDEPROBLEM = 1
  };
};

}  // namespace DirtyLittleHelpers

/**
 *  \class HierarchicalMatrix
 *  \brief Hierarchical Matrix class, which extends the EigenBase class.
 *
 *  The idea is to provide an easy to use interface to the H2-matrix
 *  from the fast boundary element method. At the moment, we inherit the
 *  traits of an Eigen::SparseMatrix, since this seems to be the minimum
 *  properties for a derived object to ensure that the matrix-vector
 *  multiplication can be specialised for HierarchicalMatrix.
 *  In particular, this allows for the use of the Eigen iterative solvers
 *  with a Hierarchical matrix.
 *
 *  \todo Maybe, we find something better then the SparsMatrix traits
 *        in the future
 **/
namespace Eigen {
/// forward definition of the HierarchicalMatrix Class in order to define traits
template <typename _Scalar>
class HierarchicalMatrix;
/// inherit the traits from the Eigen::SparseMatrix class
namespace internal {
template <typename Derived>
struct traits<HierarchicalMatrix<Derived> >
    : public internal::traits<
          SparseMatrix<typename Bembel::PDEproblemTraits<Derived>::Scalar> > {};
}  // namespace internal

// actual definition of the class
template <typename Derived>
class HierarchicalMatrix : public EigenBase<HierarchicalMatrix<Derived> > {
 public:
  // Required typedefs, constants and so on.
  /// \todo Probably, we have to change stuff here in the future
  typedef typename Bembel::PDEproblemTraits<Derived>::Scalar Scalar;
  typedef typename NumTraits<
      typename Bembel::PDEproblemTraits<Derived>::Scalar>::Real RealScalar;
  typedef Index StorageIndex;
  enum {
    ColsAtCompileTime = Dynamic,
    MaxColsAtCompileTime = Dynamic,
    IsRowMajor = false,
    Flags = NestByRefBit
  };
  // Minimum specialisation of EigenBase methods
  Index rows() const { return _H->disc->real_na; }
  Index cols() const { return _H->disc->real_na; }
  // Definition of the matrix multiplication
  template <typename Rhs>
  Product<HierarchicalMatrix, Rhs, AliasFreeProduct> operator*(
      const MatrixBase<Rhs>& x) const {
    return Product<HierarchicalMatrix, Rhs, AliasFreeProduct>(*this,
                                                              x.derived());
  }
  // Custom API:
  HierarchicalMatrix() {}
  HierarchicalMatrix(Bembel::Discretization<Derived>& disc, int np_max,
                     double eta = 1.6) {
    init_HierarchicalMatrix(disc, np_max, eta);
  }
  ~HierarchicalMatrix() {
    for (auto i = 0; i < _nct; ++i) free_cluster_tree(&_H[i]);
    if (_nct) free(_H);
  }
  void init_HierarchicalMatrix(Bembel::Discretization<Derived>& disc,
                               int np_max, double eta) {
    disc.get_pde().update_eta_by_wavenumber(eta);
    _nct = disc.get_pde()._pde.nct;
    _H = (Bembel::ct_root*)calloc(_nct, sizeof(Bembel::ct_root));
    _hmatset = get_hmatrixsettings(np_max, &(disc.get_disc()));
    _hmatset.eta = eta;
    _hmatfac.disc = &(disc.get_disc());
    _hmatfac.hmatset = &_hmatset;
    _hmatfac.assemfmats = 0;
    _hmatfac.assemsfmats = 0;
    _hmatfac.assemrkmats = 0;
    init_cluster_tree(&_hmatfac, &_H);
    _matVecHandle = Derived();
  }
  Bembel::ct_root* get_H() const { return _H; }
  const Derived& get_matVecHandle() const { return _matVecHandle; }

 private:
  int _nct;
  Bembel::ct_root* _H;
  Bembel::hmatrixsettings _hmatset; /* settings for H-matrix discretization */
  Bembel::hmatrixfactory
      _hmatfac; /* all necessary data for the matrix assembly */
  Derived _matVecHandle;
};

// Implementation of HierarchicalMatrix * Eigen::DenseVector through a
// specialization of internal::generic_product_impl:
namespace internal {
template <typename Rhs, typename Derived>
struct generic_product_impl<HierarchicalMatrix<Derived>, Rhs, SparseShape,
                            DenseShape,
                            GemvProduct>  // GEMV stands for matrix-vector
    : generic_product_impl_base<
          HierarchicalMatrix<Derived>, Rhs,
          generic_product_impl<HierarchicalMatrix<Derived>, Rhs> > {
  typedef typename Product<HierarchicalMatrix<Derived>, Rhs>::Scalar Scalar;
  template <typename Dest>
  static void scaleAndAddTo(Dest& dst, const HierarchicalMatrix<Derived>& lhs,
                            const Rhs& rhs, const Scalar& alpha) {
    // This method should implement "dst += alpha * lhs * rhs" inplace,
    // however, for iterative solvers, alpha is always equal to 1, so let's
    // not bother about it.
    assert(alpha == Scalar(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);
    // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
    // but let's do something fancier (and less efficient):
    dst += alpha * lhs.get_matVecHandle().MatVec(lhs.get_H(), rhs);
  }
};
}  // namespace internal
}  // namespace Eigen

#endif
