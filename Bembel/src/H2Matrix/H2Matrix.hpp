// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_H2MATRIX_H2MATRIX_H__
#define __BEMBEL_H2MATRIX_H2MATRIX_H__

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
  enum { YOU_CHOSE_INCOMPATIBLE_TYPES_FOR_H2MATRIX_AND_PDEPROBLEM = 1 };
};

}  // namespace DirtyLittleHelpers

/**
 *  \class H2Matrix
 *  \brief Hierarchical Matrix class, which extends the EigenBase class.
 *
 *  The idea is to provide an easy to use interface to the H2-matrix
 *  from the fast boundary element method. At the moment, we inherit the
 *  traits of an Eigen::SparseMatrix, since this seems to be the minimum
 *  properties for a derived object to ensure that the matrix-vector
 *  multiplication can be specialised for H2Matrix.
 *  In particular, this allows for the use of the Eigen iterative solvers
 *  with a Hierarchical matrix.
 *
 *  \todo Maybe, we find something better then the SparsMatrix traits
 *        in the future
 **/
namespace Eigen {
/// forward definition of the H2Matrix Class in order to define traits
template <typename ScalarT>
class H2Matrix;
/// inherit the traits from the Eigen::SparseMatrix class
namespace internal {
template <typename ScalarT>
struct traits<H2Matrix<ScalarT>>
    : public internal::traits<SparseMatrix<ScalarT>> {};
}  // namespace internal

// actual definition of the class
template <typename ScalarT>
class H2Matrix : public EigenBase<H2Matrix<ScalarT>> {
 public:
  //////////////////////////////////////////////////////////////////////////////
  /// Eigen related things
  //////////////////////////////////////////////////////////////////////////////
  // Required typedefs, constants and so on.
  typedef ScalarT Scalar;
  typedef typename NumTraits<ScalarT>::Real RealScalar;
  typedef Index StorageIndex;
  enum {
    ColsAtCompileTime = Dynamic,
    MaxColsAtCompileTime = Dynamic,
    IsRowMajor = false,
    Flags = NestByRefBit
  };
  // Minimum specialisation of EigenBase methods
  Index rows() const { return transformation_matrix_.cols(); }
  Index cols() const { return transformation_matrix_.cols(); }
  // Definition of the matrix multiplication
  template <typename Rhs>
  Product<H2Matrix, Rhs, AliasFreeProduct> operator*(
      const MatrixBase<Rhs>& x) const {
    return Product<H2Matrix, Rhs, AliasFreeProduct>(*this, x.derived());
  }
  //////////////////////////////////////////////////////////////////////////////
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  H2Matrix() {}
  /**
   * \brief Assemble H2-Matrix for linear operator linOp and AnsatzSpace
   * ansatz_space with number_of_points interpolation points in one direction of
   * the unit square (standard is number_of_points=16)
   */
  template <typename Derived>
  void init_H2Matrix(const Derived& linOp,
                     const Bembel::AnsatzSpace<Derived>& ansatz_space,
                     int number_of_points = 9) {
    // get transformation matrix from ansatz space
    transformation_matrix_ = ansatz_space.get_transformation_matrix();
    /**
     * \todo @Juergen Discuss with Michael where to initialize the parameters:
     * min_cluster_level depends on number_of_points, but this is not
     * implemented yet, also, what do the parameter mean?
     */
    // initialize block-cluster-trees and get parameters
    const int vector_dimension = Bembel::getFunctionSpaceVectorDimension<
        Bembel::LinearOperatorTraits<Derived>::Form>();
    block_cluster_tree_.resize(vector_dimension, vector_dimension);
    {
      Bembel::BlockClusterTree<Scalar> bt(linOp, ansatz_space);
      for (int i = 0; i < vector_dimension; ++i)
        for (int j = 0; j < vector_dimension; ++j)
          block_cluster_tree_(i, j) =
              bt;  //.init_BlockClusterTree(linOp, ansatz_space);
    }
    auto parameters = block_cluster_tree_(0, 0).get_parameters();
    // compute transfer and moment matrices for fmm
    int cluster_level = parameters.max_level_ - parameters.min_cluster_level_;
    if (cluster_level < 0) cluster_level = 0;
    int cluster_refinement = parameters.min_cluster_level_;
    if (cluster_refinement > parameters.max_level_)
      cluster_refinement = parameters.max_level_;
    fmm_transfer_matrices_ = Bembel::H2Multipole::computeTransferMatrices<
        Bembel::H2Multipole::ChebychevRoots>(number_of_points);
    fmm_moment_matrix_ = Bembel::H2Multipole::
        Moment2D<Bembel::H2Multipole::ChebychevRoots, Derived>::compute2DMoment(
            ansatz_space.get_superspace(), cluster_level, cluster_refinement,
            number_of_points);
    // compute interpolation points
    Eigen::VectorXd interpolation_points1D =
        Bembel::H2Multipole::ChebychevRoots(number_of_points).points_;
    Eigen::MatrixXd interpolation_points2D =
        Bembel::H2Multipole::interpolationPoints2D(interpolation_points1D);
    // compute content of tree leafs
    int polynomial_degree = ansatz_space.get_polynomial_degree();
    int polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);
    Bembel::GaussSquare<Bembel::Constants::maximum_quadrature_degree> GS;
    auto super_space = ansatz_space.get_superspace();
    auto ffield_deg = linOp.get_FarfieldQuadratureDegree(polynomial_degree);
    auto ffield_qnodes =
        Bembel::DuffyTrick::computeFfieldQnodes(super_space, GS[ffield_deg]);
    const int NumberOfFMMComponents =
        Bembel::LinearOperatorTraits<Derived>::NumberOfFMMComponents;
#pragma omp parallel
    {
#pragma omp single
      {
        // build vector of iterators
        std::vector<
            typename std::vector<Bembel::BlockClusterTree<Scalar>*>::iterator>
            leafs;
        leafs.resize(vector_dimension * vector_dimension);
        for (int i = 0; i < vector_dimension; ++i)
          for (int j = 0; j < vector_dimension; ++j)
            leafs[i * vector_dimension + j] =
                block_cluster_tree_(j, i).lbegin();
        // iterate over leafs
        for (; leafs[0] != block_cluster_tree_(0, 0).lend();) {
#pragma omp task firstprivate(leafs)
          {
            switch ((*(leafs[0]))->get_cc()) {
              // assemble dense matrix blocks
              case Bembel::BlockClusterAdmissibility::Dense: {
                const Bembel::ElementTreeNode* cluster1 =
                    (*(leafs[0]))->get_cluster1();
                const Bembel::ElementTreeNode* cluster2 =
                    (*(leafs[0]))->get_cluster2();
                const std::shared_ptr<Bembel::ElementTreeMemory> memory =
                    cluster1->get_memory();
                int block_size = std::distance(memory->cluster_begin(*cluster2),
                                               memory->cluster_end(*cluster2)) *
                                 polynomial_degree_plus_one_squared;
                std::vector<
                    Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic>>
                    F;
                for (int i = 0; i < vector_dimension * vector_dimension; ++i)
                  F.push_back(
                      Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic>(
                          block_size, block_size));
                // iterate over elements in dense matrix block
                for (auto element1 = memory->cluster_begin(*cluster1);
                     element1 != memory->cluster_end(*cluster1); ++element1) {
                  for (auto element2 = memory->cluster_begin(*cluster2);
                       element2 != memory->cluster_end(*cluster2); ++element2) {
                    Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic>
                        intval(vector_dimension *
                                   polynomial_degree_plus_one_squared,
                               vector_dimension *
                                   polynomial_degree_plus_one_squared);
                    // do integration
                    Bembel::DuffyTrick::evaluateBilinearForm(
                        linOp, super_space, *element1, *element2, GS,
                        ffield_qnodes, &intval);
                    // insert into dense matrices of all block cluster trees
                    for (int i = 0; i < vector_dimension; ++i)
                      for (int j = 0; j < vector_dimension; ++j)
                        F[i * vector_dimension + j].block(
                            polynomial_degree_plus_one_squared *
                                std::distance(memory->cluster_begin(*cluster1),
                                              element1),

                            polynomial_degree_plus_one_squared *
                                std::distance(memory->cluster_begin(*cluster2),
                                              element2),
                            polynomial_degree_plus_one_squared,
                            polynomial_degree_plus_one_squared) =
                            intval.block(j * polynomial_degree_plus_one_squared,
                                         i * polynomial_degree_plus_one_squared,
                                         polynomial_degree_plus_one_squared,
                                         polynomial_degree_plus_one_squared);
                  }
                }
                for (int i = 0; i < vector_dimension * vector_dimension; ++i)
                  (*(leafs[i]))->get_leaf().set_F(F[i]);
              } break;
              // interpolation for low-rank blocks
              case Bembel::BlockClusterAdmissibility::LowRank: {
                auto F = Bembel::H2Multipole::interpolateKernel<Derived>(
                    linOp, super_space, interpolation_points2D,
                    *((*(leafs[0]))->get_cluster1()),
                    *((*(leafs[0]))->get_cluster2()));
                for (int i = 0; i < vector_dimension; ++i) {
                  for (int j = 0; j < vector_dimension; ++j) {
                    (*(leafs[i * vector_dimension + j]))
                        ->get_leaf()
                        .set_F(F.block(j * interpolation_points2D.rows() *
                                           NumberOfFMMComponents,
                                       i * interpolation_points2D.rows() *
                                           NumberOfFMMComponents,
                                       interpolation_points2D.rows() *
                                           NumberOfFMMComponents,
                                       interpolation_points2D.rows() *
                                           NumberOfFMMComponents));
                    (*(leafs[i * vector_dimension + j]))
                        ->get_leaf()
                        .set_low_rank_flag(true);
                  }
                }
              } break;
              // this leaf is not a low-rank block and not a dense block,
              // thus it is not a leaf -> error
              default:
                assert(0 && "This should never happen");
                break;
            }
          }
          // increment leafs of all block-cluster trees
          for (auto leafsit = leafs.begin(); leafsit != leafs.end(); ++leafsit)
            ++(*leafsit);
        }
      }
    }
  }
  //////////////////////////////////////////////////////////////////////////////
  /// methods
  //////////////////////////////////////////////////////////////////////////////
  Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> get_dense() const {
    Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> dense(rows(),
                                                                 cols());
    for (int i = 0; i < cols(); ++i) {
      Eigen::VectorXd unit(cols());
      unit.setZero();
      unit(i) = 1.;
      dense.col(i) = (*this) * unit;
    }
    return dense;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// getter
  //////////////////////////////////////////////////////////////////////////////
  const Eigen::SparseMatrix<double> get_transformation_matrix() const {
    return transformation_matrix_;
  }
  const Eigen::MatrixXd get_fmm_transfer_matrices() const {
    return fmm_transfer_matrices_;
  }
  const std::vector<Eigen::MatrixXd> get_fmm_moment_matrix() const {
    return fmm_moment_matrix_;
  }
  const Bembel::GenericMatrix<Bembel::BlockClusterTree<ScalarT>>&
  get_block_cluster_tree() const {
    return block_cluster_tree_;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// private members
  //////////////////////////////////////////////////////////////////////////////
 private:
  // we declare functionality which has not been implemented (yet)
  // to be private
  H2Matrix(const H2Matrix<ScalarT>& H);
  H2Matrix(H2Matrix<ScalarT>&& H);
  H2Matrix& operator=(const H2Matrix<ScalarT>& H);
  H2Matrix& operator=(H2Matrix<ScalarT>&& H);
  Eigen::SparseMatrix<double> transformation_matrix_;
  Bembel::GenericMatrix<Bembel::BlockClusterTree<ScalarT>> block_cluster_tree_;
  Eigen::MatrixXd fmm_transfer_matrices_;
  std::vector<Eigen::MatrixXd> fmm_moment_matrix_;
};  // namespace Eigen

/**
 * \brief Implementation of H2Matrix * Eigen::DenseVector through a
 * specialization of internal::generic_product_impl
 */
namespace internal {
template <typename Rhs, typename ScalarT>
struct generic_product_impl<H2Matrix<ScalarT>, Rhs, SparseShape, DenseShape,
                            GemvProduct>  // GEMV stands for matrix-vector
    : generic_product_impl_base<H2Matrix<ScalarT>, Rhs,
                                generic_product_impl<H2Matrix<ScalarT>, Rhs>> {
  template <typename Dest>
  static void scaleAndAddTo(Dest& dst, const H2Matrix<ScalarT>& lhs,
                            const Rhs& rhs, const ScalarT& alpha) {
    // This method should implement "dst += alpha * lhs * rhs" inplace, however,
    // for iterative solvers, alpha is always equal to 1, so let's not bother
    // about it.
    assert(alpha == ScalarT(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);

    // get H2-data
    int max_level =
        lhs.get_block_cluster_tree()(0, 0).get_parameters().max_level_;
    int min_cluster_level =
        lhs.get_block_cluster_tree()(0, 0).get_parameters().min_cluster_level_;
    auto moment_matrix = lhs.get_fmm_moment_matrix();
    auto transfer_matrices = lhs.get_fmm_transfer_matrices();
    int vector_dimension = moment_matrix.size();

    // go discontinuous in rhs
    Matrix<ScalarT, Dynamic, 1> long_rhs_all =
        (lhs.get_transformation_matrix() * rhs).eval();
    int vector_component_size = long_rhs_all.rows() / vector_dimension;

    // initialize destination
    Matrix<ScalarT, Dynamic, 1> long_dst_all(long_rhs_all.rows());
    long_dst_all.setZero();

    for (int col_component = 0; col_component < vector_dimension;
         ++col_component) {
      for (int row_component = 0; row_component < vector_dimension;
           ++row_component) {
        Matrix<ScalarT, Dynamic, 1> long_rhs = long_rhs_all.segment(
            col_component * vector_component_size, vector_component_size);
        Matrix<ScalarT, Dynamic, 1> long_dst(long_rhs.rows());
        long_dst.setZero();

        // split long rhs into pieces by reshaping
        Matrix<ScalarT, Dynamic, Dynamic> long_rhs_matrix =
            Map<Matrix<ScalarT, Dynamic, Dynamic>>(
                long_rhs.data(), moment_matrix[col_component].cols(),
                long_rhs.rows() / moment_matrix[col_component].cols());

        // do forward-transformation
        std::vector<Matrix<ScalarT, Dynamic, Dynamic>> long_rhs_forward =
            Bembel::H2Multipole::forwardTransformation(
                moment_matrix[col_component], transfer_matrices,
                max_level - min_cluster_level, long_rhs_matrix);

#pragma omp parallel
        {
          // initialize target for each process
          Matrix<ScalarT, Dynamic, 1> my_long_dst(long_dst.rows());

          // initialize target of backward-transformation
          std::vector<Matrix<ScalarT, Dynamic, Dynamic>> my_long_dst_backward;
          for (int i = 0; i < long_rhs_forward.size(); ++i)
            my_long_dst_backward.push_back(
                Matrix<ScalarT, Dynamic, Dynamic>::Zero(
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
                  const std::shared_ptr<Bembel::ElementTreeMemory> memory =
                      cluster1->get_memory();
                  int cluster1_col = cluster1->id_;
                  int cluster2_col = cluster2->id_;
                  my_long_dst_backward[fmm_level].col(cluster1_col) +=
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
    dst += lhs.get_transformation_matrix().transpose() * long_dst_all;
  }
};

}  // namespace internal
}  // namespace Eigen

#endif
