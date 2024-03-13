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
#ifndef BEMBEL_SRC_H2MATRIX_H2MATRIX_HPP_
#define BEMBEL_SRC_H2MATRIX_H2MATRIX_HPP_

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
 **/
namespace Eigen {
/// forward definition of the H2Matrix Class in order to define traits
template <typename ScalarT>
class H2Matrix;
/// small modification of internal traits of SparseMatrix
namespace internal {
template <typename ScalarT>
struct traits<H2Matrix<ScalarT>> {
  typedef ScalarT Scalar;
  typedef int StorageIndex;
  typedef H2 StorageKind;
  typedef MatrixXpr XprKind;
  enum {
    RowsAtCompileTime = Dynamic,
    ColsAtCompileTime = Dynamic,
    MaxRowsAtCompileTime = Dynamic,
    MaxColsAtCompileTime = Dynamic,
    Flags = NestByRefBit
  };
};

// this struct is necessary for compatibility with iterative solvers
template <typename ScalarT>
struct is_ref_compatible<H2Matrix<ScalarT>> {
  enum { value = false };
};
template <typename ScalarT>
struct is_ref_compatible<const H2Matrix<ScalarT>>
    : is_ref_compatible<H2Matrix<ScalarT>> {};

}  // namespace internal

// actual definition of the class
template <typename ScalarT>
class H2Matrix : public H2MatrixBase<H2Matrix<ScalarT>> {
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
   * the unit square (standard is number_of_points=9)
   */
  template <typename Derived>
  void init_H2Matrix(const Derived& linOp,
                     const Bembel::AnsatzSpace<Derived>& ansatz_space,
                     int number_of_points = 9) {
    // get transformation matrix from ansatz space
    transformation_matrix_ = ansatz_space.get_transformation_matrix();
    /**
     * \todo Juergen Discuss with Michael where to initialize the parameters:
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
              bt;  // .init_BlockClusterTree(linOp, ansatz_space);
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
    std::vector<ElementSurfacePoints> ffield_qnodes =
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
                int block_size =
                    std::distance(cluster2->begin(), cluster2->end()) *
                    polynomial_degree_plus_one_squared;
                std::vector<
                    Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic>>
                    F;
                for (int i = 0; i < vector_dimension * vector_dimension; ++i)
                  F.push_back(
                      Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic>(
                          block_size, block_size));
                // iterate over elements in dense matrix block
                unsigned int cl1index = 0;
                unsigned int cl2index = 0;
                for (const auto& element1 : *cluster1) {
                  cl2index = 0;
                  for (const auto& element2 : *cluster2) {
                    Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic>
                        intval(vector_dimension *
                                   polynomial_degree_plus_one_squared,
                               vector_dimension *
                                   polynomial_degree_plus_one_squared);
                    // do integration
                    Bembel::DuffyTrick::evaluateBilinearForm(
                        linOp, super_space, element1, element2, GS,
                        ffield_qnodes[element1.id_],
                        ffield_qnodes[element2.id_], &intval);
                    // insert into dense matrices of all block cluster trees
                    for (int i = 0; i < vector_dimension; ++i)
                      for (int j = 0; j < vector_dimension; ++j)
                        F[i * vector_dimension + j].block(
                            polynomial_degree_plus_one_squared * cl1index,
                            polynomial_degree_plus_one_squared * cl2index,
                            polynomial_degree_plus_one_squared,
                            polynomial_degree_plus_one_squared) =
                            intval.block(j * polynomial_degree_plus_one_squared,
                                         i * polynomial_degree_plus_one_squared,
                                         polynomial_degree_plus_one_squared,
                                         polynomial_degree_plus_one_squared);
                    ++cl2index;
                  }
                  ++cl1index;
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
};

namespace internal {

// adaption from SparseMatrix from Eigen
template <typename Scalar>
struct evaluator<H2Matrix<Scalar>> : evaluator<H2MatrixBase<H2Matrix<Scalar>>> {
  typedef evaluator<H2MatrixBase<H2Matrix<Scalar>>> Base;
  typedef H2Matrix<Scalar> H2MatrixType;
  evaluator() : Base() {}
  explicit evaluator(const H2MatrixType& mat) : Base(mat) {}
};

}  // namespace internal
}  // namespace Eigen

#endif  // BEMBEL_SRC_H2MATRIX_H2MATRIX_HPP_
