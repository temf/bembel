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
#ifndef BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2DENSEPRODUCT_HPP_
#define BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2DENSEPRODUCT_HPP_

// The contents of this file are a modification of the file
// SparseDenseProduct.h from the Eigen library
namespace Eigen {

namespace internal {

template <typename H2LhsType, typename DenseRhsType, typename DenseResType,
          typename AlphaType,
          int LhsStorageOrder =
              ((H2LhsType::Flags & RowMajorBit) == RowMajorBit) ? RowMajor
                                                                : ColMajor,
          bool ColPerCol = ((DenseRhsType::Flags & RowMajorBit) == 0) ||
                           DenseRhsType::ColsAtCompileTime == 1>
struct H2_time_dense_product_impl;

/**
 * \brief H2-matrix-vector multiplication, works also for matrices by iterating
 * over the columns.
 */
template <typename ScalarH2, typename DenseRhsType, typename DenseResType,
          typename AlphaType>
struct H2_time_dense_product_impl<H2Matrix<ScalarH2>, DenseRhsType,
                                  DenseResType, AlphaType, ColMajor, true> {
  typedef typename traits<DenseRhsType>::Scalar ScalarRhs;
  typedef typename traits<DenseResType>::Scalar ScalarRes;
  static void run(const H2Matrix<ScalarH2>& lhs, const DenseRhsType& rhs,
                  DenseResType& res, const AlphaType& alpha) {
    // get H2-data
    int max_level =
        lhs.get_block_cluster_tree()(0, 0).get_parameters().max_level_;
    int min_cluster_level =
        lhs.get_block_cluster_tree()(0, 0).get_parameters().min_cluster_level_;
    auto moment_matrix = lhs.get_fmm_moment_matrix();
    auto transfer_matrices = lhs.get_fmm_transfer_matrices();
    int vector_dimension = moment_matrix.size();

    for (Index c = 0; c < rhs.cols(); ++c) {
      // go discontinuous in rhs
      Matrix<ScalarH2, Dynamic, 1> long_rhs_all =
          (lhs.get_transformation_matrix() * rhs.col(c)).eval();
      int vector_component_size = long_rhs_all.rows() / vector_dimension;

      // initialize destination
      Matrix<ScalarRes, Dynamic, 1> long_dst_all(long_rhs_all.rows());
      long_dst_all.setZero();

      for (int col_component = 0; col_component < vector_dimension;
           ++col_component) {
        for (int row_component = 0; row_component < vector_dimension;
             ++row_component) {
          Matrix<ScalarRhs, Dynamic, 1> long_rhs = long_rhs_all.segment(
              col_component * vector_component_size, vector_component_size);
          Matrix<ScalarRhs, Dynamic, 1> long_dst(long_rhs.rows());
          long_dst.setZero();

          // split long rhs into pieces by reshaping
          Matrix<ScalarRhs, Dynamic, Dynamic> long_rhs_matrix =
              Map<Matrix<ScalarRhs, Dynamic, Dynamic>>(
                  long_rhs.data(), moment_matrix[col_component].cols(),
                  long_rhs.rows() / moment_matrix[col_component].cols());

          // do forward-transformation
          std::vector<Matrix<ScalarRhs, Dynamic, Dynamic>> long_rhs_forward =
              Bembel::H2Multipole::forwardTransformation(
                  moment_matrix[col_component], transfer_matrices,
                  max_level - min_cluster_level, long_rhs_matrix);

#pragma omp parallel
          {
            // initialize target for each process
            Matrix<ScalarRes, Dynamic, 1> my_long_dst(long_dst.rows());

            // initialize target of backward-transformation
            std::vector<Matrix<ScalarRes, Dynamic, Dynamic>>
                my_long_dst_backward;
            for (int i = 0; i < long_rhs_forward.size(); ++i)
              my_long_dst_backward.push_back(
                  Matrix<ScalarRes, Dynamic, Dynamic>::Zero(
                      long_rhs_forward[i].rows(), long_rhs_forward[i].cols()));
            my_long_dst.setZero();

            // matrix-vector
            for (auto leaf =
                     lhs.get_block_cluster_tree()(row_component, col_component)
                         .clbegin();
                 leaf !=
                 lhs.get_block_cluster_tree()(row_component, col_component)
                     .clend();
                 ++leaf) {
#pragma omp single nowait
              {
                switch ((*leaf)->get_cc()) {
                  // deal with matrix blocks
                  case Bembel::BlockClusterAdmissibility::Dense: {
                    my_long_dst.segment((*leaf)->get_row_start_index(),
                                        (*leaf)->get_row_end_index() -
                                            (*leaf)->get_row_start_index()) +=
                        (*leaf)->get_leaf().get_F() *
                        long_rhs.segment((*leaf)->get_col_start_index(),
                                         (*leaf)->get_col_end_index() -
                                             (*leaf)->get_col_start_index());
                  } break;
                  // deal with low-rank blocks
                  case Bembel::BlockClusterAdmissibility::LowRank: {
                    const Bembel::ElementTreeNode* cluster1 =
                        (*leaf)->get_cluster1();
                    const Bembel::ElementTreeNode* cluster2 =
                        (*leaf)->get_cluster2();
                    int cluster_level = cluster1->get_level();
                    int fmm_level = lhs.get_block_cluster_tree()(0, 0)
                                        .get_parameters()
                                        .max_level_ -
                                    lhs.get_block_cluster_tree()(0, 0)
                                        .get_parameters()
                                        .min_cluster_level_ -
                                    (*leaf)->get_cluster1()->get_level();
                    int cluster1_col = cluster1->id_;
                    int cluster2_col = cluster2->id_;
                    my_long_dst_backward[fmm_level].col(cluster1_col) +=
                        (*leaf)->get_leaf().get_F() *
                        long_rhs_forward[fmm_level].col(cluster2_col);
                  } break;
                  // this leaf is not a low-rank block and not a dense block,
                  // thus it is not a leaf -> error
                  default:
                    assert(0 && "This should never happen");
                    break;
                }
              }
            }

            // do backward transformation
            my_long_dst += Bembel::H2Multipole::backwardTransformation(
                moment_matrix[row_component], transfer_matrices,
                max_level - min_cluster_level, my_long_dst_backward);

#pragma omp critical
            long_dst += my_long_dst;
          }

          // finish off vector component
          long_dst_all.segment(row_component * vector_component_size,
                               vector_component_size) += long_dst;
        }
      }

      // go continuous and write output
      res.col(c) +=
          alpha * lhs.get_transformation_matrix().transpose() * long_dst_all;
    }
  }
};

template <typename H2LhsType, typename DenseRhsType, typename DenseResType,
          typename AlphaType>
inline void H2_time_dense_product(const H2LhsType& lhs, const DenseRhsType& rhs,
                                  DenseResType& res, const AlphaType& alpha) {
  H2_time_dense_product_impl<H2LhsType, DenseRhsType, DenseResType,
                             AlphaType>::run(lhs, rhs, res, alpha);
}

/**
 * \brief overwrite eigen implementation by using distributive law, i.e.,
 *compute A*c+B*c instead of (A+B)*c
 */
template <typename BinaryOp, typename BinaryLhs, typename BinaryRhs,
          typename Rhs, int ProductType>
struct generic_product_impl<CwiseBinaryOp<BinaryOp, BinaryLhs, BinaryRhs>, Rhs,
                            H2, DenseShape, ProductType>
    : generic_product_impl_base<
          CwiseBinaryOp<BinaryOp, BinaryLhs, BinaryRhs>, Rhs,
          generic_product_impl<CwiseBinaryOp<BinaryOp, BinaryLhs, BinaryRhs>,
                               Rhs, H2, DenseShape, ProductType>> {
  typedef CwiseBinaryOp<BinaryOp, BinaryLhs, BinaryRhs> Lhs;
  typedef typename Product<Lhs, Rhs>::Scalar Scalar;

  // we only need to specify scaleAndAddTo, everything else is build on this in
  // the background
  template <typename Dest>
  static void scaleAndAddTo(Dest& dst, const Lhs& lhs, const Rhs& rhs,
                            const Scalar& alpha) {
    typedef Product<BinaryLhs, Rhs, ProductType> LeftProduct;
    typedef Product<BinaryRhs, Rhs, ProductType> RightProduct;
    typedef CwiseBinaryOp<BinaryOp, LeftProduct, RightProduct> NewCwiseBinaryOp;
    typedef Matrix<Scalar, Dynamic, Dynamic> DenseMatrix;
    typedef CwiseNullaryOp<scalar_constant_op<Scalar>, DenseMatrix>
        ScalarCwiseNullaryOp;
    typedef scalar_constant_op<Scalar> ScalarConstantOp;
    typedef scalar_product_op<Scalar, Scalar> ScalarTimesOp;
    typedef CwiseBinaryOp<ScalarTimesOp, ScalarCwiseNullaryOp, NewCwiseBinaryOp>
        ScalarTimesNewCwiseBinaryOp;
    typedef typename traits<Lhs>::Scalar LhsScalar;
    typedef typename traits<Rhs>::Scalar RhsScalar;
    typedef add_assign_op<LhsScalar, RhsScalar> addAssignOp;

    static_assert(
        is_same<BinaryOp, scalar_sum_op<LhsScalar, RhsScalar>>::value ||
            is_same<BinaryOp,
                    scalar_difference_op<LhsScalar, RhsScalar>>::value,
        "Product of CwiseBinaryOp not defined for this BinaryOp");

    LeftProduct lprod(lhs.lhs(), rhs);
    RightProduct rprod(lhs.rhs(), rhs);
    NewCwiseBinaryOp xpr_prod(lprod, rprod, lhs.functor());
    ScalarCwiseNullaryOp xpr_scalar(lprod.rows(), lprod.cols(),
                                    ScalarConstantOp(alpha));
    ScalarTimesNewCwiseBinaryOp xpr(xpr_scalar, xpr_prod, ScalarTimesOp());
    Assignment<Dest, ScalarTimesNewCwiseBinaryOp, addAssignOp,
               Dense2Dense>::run(dst, xpr, addAssignOp());
  }
};
/**
 * \brief overwrite eigen implementation by using associative law to compute
 *a*(H*x) instead of (a*H)*x
 */
template <typename ProdLhsScalar, typename ProdRhsScalar, typename BinaryLhs,
          typename BinaryRhs, typename Rhs, int ProductType>
struct generic_product_impl<
    CwiseBinaryOp<scalar_product_op<ProdLhsScalar, ProdRhsScalar>, BinaryLhs,
                  BinaryRhs>,
    Rhs, H2, DenseShape, ProductType>
    : generic_product_impl_base<
          CwiseBinaryOp<scalar_product_op<ProdLhsScalar, ProdRhsScalar>,
                        BinaryLhs, BinaryRhs>,
          Rhs,
          generic_product_impl<
              CwiseBinaryOp<scalar_product_op<ProdLhsScalar, ProdRhsScalar>,
                            BinaryLhs, BinaryRhs>,
              Rhs, H2, DenseShape, ProductType>> {
  typedef CwiseBinaryOp<scalar_product_op<ProdLhsScalar, ProdRhsScalar>,
                        BinaryLhs, BinaryRhs>
      Lhs;
  typedef typename Product<Lhs, Rhs>::Scalar Scalar;

  // we only need to specify scaleAndAddTo, everything else is build on this in
  // the background
  template <typename Dest>
  static void scaleAndAddTo(Dest& dst, const Lhs& lhs, const Rhs& rhs,
                            const Scalar& alpha) {
    Scalar beta = lhs.lhs().functor()();
    generic_product_impl<BinaryRhs, Rhs, H2, DenseShape,
                         ProductType>::scaleAndAddTo(dst, lhs.rhs(), rhs,
                                                     alpha * beta);
  }
};
// same overwrite, but for const CwiseBinaryOp by inheritance from non-const
// version
template <typename BinaryOp, typename BinaryLhs, typename BinaryRhs,
          typename Rhs, int ProductType>
struct generic_product_impl<const CwiseBinaryOp<BinaryOp, BinaryLhs, BinaryRhs>,
                            Rhs, H2, DenseShape, ProductType>
    : generic_product_impl<CwiseBinaryOp<BinaryOp, BinaryLhs, BinaryRhs>, Rhs,
                           H2, DenseShape, ProductType> {};

template <typename Lhs, typename Rhs, int ProductType>
struct generic_product_impl<Lhs, Rhs, H2, DenseShape, ProductType>
    : generic_product_impl_base<
          Lhs, Rhs,
          generic_product_impl<Lhs, Rhs, H2, DenseShape, ProductType>> {
  typedef typename Product<Lhs, Rhs>::Scalar Scalar;

  // we only need to specify scaleAndAddTo, everything else is build on this in
  // the background
  template <typename Dest>
  static void scaleAndAddTo(Dest& dst, const Lhs& lhs, const Rhs& rhs,
                            const Scalar& alpha) {
    typedef
        typename nested_eval<Lhs, ((Rhs::Flags & RowMajorBit) == 0)
                                      ? 1
                                      : Rhs::ColsAtCompileTime>::type LhsNested;
    typedef typename nested_eval<
        Rhs, ((Lhs::Flags & RowMajorBit) == 0) ? 1 : Dynamic>::type RhsNested;
    LhsNested lhsNested(lhs);
    RhsNested rhsNested(rhs);
    internal::H2_time_dense_product(lhsNested, rhsNested, dst, alpha);
  }
};

template <typename Lhs, typename Rhs, int Options, int ProductTag>
struct product_evaluator<Product<Lhs, Rhs, Options>, ProductTag, H2, DenseShape>
    : public evaluator<typename Product<Lhs, Rhs, Options>::PlainObject> {
  typedef Product<Lhs, Rhs, Options> XprType;
  typedef typename XprType::PlainObject PlainObject;
  typedef evaluator<PlainObject> Base;

  enum { Flags = Base::Flags };

  explicit product_evaluator(const XprType& xpr)
      : m_result(xpr.rows(), xpr.cols()) {
    ::new (static_cast<Base*>(this)) Base(m_result);
    generic_product_impl<Lhs, Rhs, H2, DenseShape, ProductTag>::evalTo(
        m_result, xpr.lhs(), xpr.rhs());
  }

 protected:
  PlainObject m_result;
};

}  // end namespace internal

}  // end namespace Eigen

#endif  // BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2DENSEPRODUCT_HPP_
