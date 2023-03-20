// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2DENSEPRODUCT_H__
#define __BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2DENSEPRODUCT_H__

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

template <typename H2LhsType, typename DenseRhsType, typename DenseResType,
          typename AlphaType>
struct H2_time_dense_product_impl<H2LhsType, DenseRhsType, DenseResType,
                                  AlphaType, ColMajor, true> {
  static void run(const H2LhsType& lhs, const DenseRhsType& rhs,
                  DenseResType& res, const AlphaType& alpha) {
    typedef typename internal::remove_all<H2LhsType>::type Lhs;
    typedef typename internal::remove_all<DenseRhsType>::type Rhs;
    typedef typename internal::remove_all<DenseResType>::type Res;
    typedef Matrix<Lhs, Dynamic, 1> RhsVec;
    typedef Matrix<Lhs, Dynamic, 1> ResVec;
    typedef Matrix<Lhs, Dynamic, Dynamic> RhsMat;
    typedef Matrix<Lhs, Dynamic, Dynamic> ResMat;

    // get H2-data
    const int max_level =
        lhs.get_block_cluster_tree()(0, 0).get_parameters().max_level_;
    const int min_cluster_level =
        lhs.get_block_cluster_tree()(0, 0).get_parameters().min_cluster_level_;
    const std::vector<MatrixXd>& moment_matrix = lhs.get_fmm_moment_matrix();
    const MatrixXd& transfer_matrices = lhs.get_fmm_transfer_matrices();
    int vector_dimension = moment_matrix.size();

    // go discontinuous in rhs
    RhsVec long_rhs_all = (lhs.get_transformation_matrix() * rhs).eval();
    int vector_component_size = long_rhs_all.rows() / vector_dimension;

    // initialize destination
    RhsVec long_res_all(long_rhs_all.rows());
    long_res_all.setZero();

    for (int col_component = 0; col_component < vector_dimension;
         ++col_component) {
      for (int row_component = 0; row_component < vector_dimension;
           ++row_component) {
        RhsVec long_rhs = long_rhs_all.segment(
            col_component * vector_component_size, vector_component_size);
        ResVec long_res(long_rhs.rows());
        long_res.setZero();

        // split long rhs into pieces by reshaping
        RhsMat long_rhs_matrix =
            Map<RhsMat>(long_rhs.data(), moment_matrix[col_component].cols(),
                        long_rhs.rows() / moment_matrix[col_component].cols());

        // do forward-transformation
        std::vector<RhsMat> long_rhs_forward =
            Bembel::H2Multipole::forwardTransformation(
                moment_matrix[col_component], transfer_matrices,
                max_level - min_cluster_level, long_rhs_matrix);

#pragma omp parallel
        {
          // initialize target for each process
          ResVec my_long_res(long_res.rows());

          // initialize target of backward-transformation
          std::vector<ResMat> my_long_res_backward;
          for (int i = 0; i < long_rhs_forward.size(); ++i)
            my_long_res_backward.push_back(ResMat::Zero(
                long_rhs_forward[i].rows(), long_rhs_forward[i].cols()));
          my_long_res.setZero();

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
                  my_long_res.segment((*leaf)->get_row_start_index(),
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
                  my_long_res_backward[fmm_level].col(cluster1_col) +=
                      (*leaf)->get_leaf().get_F() *
                      long_rhs_forward[fmm_level].col(cluster2_col);
                } break;
                // this leaf is not a low-rank block and not a dense block, thus
                // it is not a leaf -> error
                default:
                  assert(0 && "This should never happen");
                  break;
              }
            }
          }

          // do backward transformation
          my_long_res += Bembel::H2Multipole::backwardTransformation(
              moment_matrix[row_component], transfer_matrices,
              max_level - min_cluster_level, my_long_res_backward);

#pragma omp critical
          long_res += my_long_res;
        }

        // finish off vector component
        long_res_all.segment(row_component * vector_component_size,
                             vector_component_size) += long_res;
      }
    }

    // go continuous and write output
    res += alpha * lhs.get_transformation_matrix().transpose() * long_res_all;
  }
};

template <typename H2LhsType, typename DenseRhsType, typename DenseResType,
          typename AlphaType>
inline void H2_time_dense_product(const H2LhsType& lhs, const DenseRhsType& rhs,
                                  DenseResType& res, const AlphaType& alpha) {
  H2_time_dense_product_impl<H2LhsType, DenseRhsType, DenseResType,
                             AlphaType>::run(lhs, rhs, res, alpha);
}

template <typename Lhs, typename Rhs, int ProductType>
struct generic_product_impl<Lhs, Rhs, H2, DenseShape, ProductType>
    : generic_product_impl_base<
          Lhs, Rhs,
          generic_product_impl<Lhs, Rhs, H2, DenseShape, ProductType>> {
  typedef typename Product<Lhs, Rhs>::Scalar Scalar;

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

}  // end namespace internal

}  // end namespace Eigen

#endif  // __BEMBEL_SRC_H2MATRIX_EIGENHELPER_H2ASSIGN_H__