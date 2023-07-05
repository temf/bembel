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
#ifndef BEMBEL_SRC_H2MATRIX_BLOCKCLUSTERTREE_HPP_
#define BEMBEL_SRC_H2MATRIX_BLOCKCLUSTERTREE_HPP_

namespace Bembel {

template <typename Scalar>
class BlockClusterTree;

struct BlockClusterAdmissibility {
  enum { Refine = 0, LowRank = 1, Dense = 2 };
};

template <typename Scalar>
struct BlockClusterTreeParameters {
  BlockClusterTreeParameters()
      : eta_(-1), min_cluster_level_(-1), max_level_(-1) {}
  GaussSquare<Constants::maximum_quadrature_degree> GS_;
  Eigen::Matrix<double, 12, Eigen::Dynamic> ffield_qnodes_;
  double eta_;             // eta from admissibility condition
  int ffield_deg_;         // todo @Michael comment this
  int min_cluster_level_;  // a tree leaf has 4^min_cluster_level_ elements
  int max_level_;          // depth of block cluster tree
  int polynomial_degree_;  // todo @Michael comment this
  int polynomial_degree_plus_one_squared_;  // todo @Michael comment this
  const ElementTreeNode *et_root_;
};

template <typename Scalar>
class BlockClusterTree {
 public:
  //////////////////////////////////////////////////////////////////////////////
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  BlockClusterTree()
      : parameters_(nullptr),
        sons_(GenericMatrix<BlockClusterTree>()),
        leaf_(nullptr),
        leaf_pointers_(nullptr),
        cluster1_(nullptr),
        cluster2_(nullptr),
        rows_(0),
        cols_(0),
        cc_(-1) {}

  BlockClusterTree(BlockClusterTree &&other) noexcept {
    parameters_ = other.parameters_;
    sons_ = std::move(other.sons_);
    leaf_ = other.leaf_;
    leaf_pointers_ = other.leaf_pointers_;
    cluster1_ = other.cluster1_;
    cluster2_ = other.cluster2_;
    rows_ = other.rows_;
    cols_ = other.cols_;
    cc_ = other.cc_;
    updateLeafPointers();
  }

  BlockClusterTree(const BlockClusterTree &other) noexcept {
    parameters_ = other.parameters_;
    sons_ = other.sons_;
    // tree leaf is deep copied for const copy constructor! This has to be
    // taken care of in some sense when going for Hmatrix arithmetics...
    if (other.leaf_) {
      leaf_ = std::make_shared<
          TreeLeaf<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>>();
      *leaf_ = static_cast<const TreeLeaf<
          Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> &>(
          *(other.leaf_));
    }
    leaf_pointers_ = std::make_shared<std::vector<BlockClusterTree *>>();
    cluster1_ = other.cluster1_;
    cluster2_ = other.cluster2_;
    rows_ = other.rows_;
    cols_ = other.cols_;
    cc_ = other.cc_;
    updateLeafPointers();
  }

  template <typename Derived>
  BlockClusterTree(const LinearOperatorBase<Derived> &linear_operator,
                   const AnsatzSpace<Derived> &ansatz_space) {
    init_BlockClusterTree(linear_operator, ansatz_space);
  }
  //////////////////////////////////////////////////////////////////////////////
  /// assignment operator
  //////////////////////////////////////////////////////////////////////////////
  BlockClusterTree &operator=(BlockClusterTree other) noexcept {
    std::swap(parameters_, other.parameters_);
    std::swap(sons_, other.sons_);
    std::swap(leaf_, other.leaf_);
    std::swap(leaf_pointers_, other.leaf_pointers_);
    std::swap(cluster1_, other.cluster1_);
    std::swap(cluster2_, other.cluster2_);
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(cc_, other.cc_);
    updateLeafPointers();
    return *this;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// init
  //////////////////////////////////////////////////////////////////////////////
  void set_parameters(double eta = 1.6, int min_cluster_level = 1) {
    parameters_ = std::make_shared<BlockClusterTreeParameters<Scalar>>();
    parameters_->eta_ = eta;
    parameters_->min_cluster_level_ = min_cluster_level;
    return;
  }

  template <typename Derived>
  void init_BlockClusterTree(const LinearOperatorBase<Derived> &linOp,
                             const AnsatzSpace<Derived> &ansatz_space) {
    if (parameters_ == nullptr) set_parameters();
    // get element tree from ansatz space
    const ElementTree &element_tree =
        ansatz_space.get_superspace().get_mesh().get_element_tree();
    parameters_->et_root_ = std::addressof(element_tree.root());
    cluster1_ = std::addressof(element_tree.root());
    cluster2_ = std::addressof(element_tree.root());
    // set parameters for matrix assembly
    parameters_->max_level_ = element_tree.get_max_level();
    parameters_->polynomial_degree_ =
        ansatz_space.get_superspace().get_polynomial_degree();
    parameters_->polynomial_degree_plus_one_squared_ =
        (parameters_->polynomial_degree_ + 1) *
        (parameters_->polynomial_degree_ + 1);
    parameters_->ffield_deg_ =
        linOp.get_FarfieldQuadratureDegree(parameters_->polynomial_degree_);
    // set up leaf_pointers
    leaf_pointers_ = std::make_shared<std::vector<BlockClusterTree *>>();
    leaf_pointers_->clear();
    // block cluster tree assembly
    rows_ = element_tree.get_number_of_elements();
    cols_ = element_tree.get_number_of_elements();
    // we let appendSubtree handle everything, since the root always
    // returns 0 in compare cluster
    appendSubtree(linOp, ansatz_space, cluster1_, cluster2_);
    updateLeafPointers();
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// methods
  //////////////////////////////////////////////////////////////////////////////
  template <typename Derived>
  void appendSubtree(const LinearOperatorBase<Derived> &linear_operator,
                     const AnsatzSpace<Derived> &ansatz_space,
                     const ElementTreeNode *cluster1,
                     const ElementTreeNode *cluster2) {
    cc_ = compareCluster(*cluster1, *cluster2);
    // there are children to handle
    if (cc_ == BlockClusterAdmissibility::Refine) {
      // reserve memory for sons_
      sons_.resize(cluster1->sons_.size(), cluster2->sons_.size());
      for (auto j = 0; j < sons_.cols(); ++j)
        for (auto i = 0; i < sons_.rows(); ++i) {
          const ElementTreeNode &son1 = cluster1->sons_[i];
          const ElementTreeNode &son2 = cluster2->sons_[j];
          sons_(i, j).parameters_ = parameters_;
          sons_(i, j).cluster1_ = std::addressof(son1);
          sons_(i, j).cluster2_ = std::addressof(son2);
          sons_(i, j).rows_ = std::distance(son1.begin(), son1.end());
          sons_(i, j).cols_ = std::distance(son2.begin(), son2.end());
          // let recursion handle the rest
          sons_(i, j).appendSubtree(linear_operator, ansatz_space,
                                    std::addressof(son1), std::addressof(son2));
        }
    } else {
      leaf_ = std::make_shared<
          TreeLeaf<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>>();
    }
    return;
  }

  void printNode() const {
    std::cout << "{" << std::endl;
    std::cout << "cluster: " << cluster1_ << " & " << cluster2_ << std::endl;
    std::cout << "rows: " << rows_ << " cols: " << cols_ << std::endl;
    std::cout << "format: " << cc_ << std::endl;
    std::cout << "children: " << sons_.rows() << " x " << sons_.cols()
              << std::endl;
    std::cout << "leaf: " << leaf_ << std::endl;
    std::cout << "parameters: " << parameters_ << std::endl;

    std::cout << "}" << std::endl;
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// getter
  //////////////////////////////////////////////////////////////////////////////
  int get_cc() const { return cc_; }
  int rows() const { return rows_; }
  int cols() const { return cols_; }
  const ElementTreeNode *get_cluster1() const { return cluster1_; }
  const ElementTreeNode *get_cluster2() const { return cluster2_; }
  int get_level() const { return get_cluster1()->get_level(); }
  typename std::vector<BlockClusterTree *>::iterator lbegin() {
    return leaf_pointers_->begin();
  }
  typename std::vector<BlockClusterTree *>::iterator lend() {
    return leaf_pointers_->end();
  }
  typename std::vector<BlockClusterTree *>::const_iterator clbegin() const {
    return leaf_pointers_->begin();
  }
  typename std::vector<BlockClusterTree *>::const_iterator clend() const {
    return leaf_pointers_->end();
  }
  int get_row_start_index() {
    return get_parameters().polynomial_degree_plus_one_squared_ *
           std::distance(get_parameters().et_root_->begin(),
                         cluster1_->begin());
  }
  int get_row_end_index() {
    return get_parameters().polynomial_degree_plus_one_squared_ *
           std::distance(get_parameters().et_root_->begin(), cluster1_->end());
  }
  int get_col_start_index() {
    return get_parameters().polynomial_degree_plus_one_squared_ *
           std::distance(get_parameters().et_root_->begin(),
                         cluster2_->begin());
  }
  int get_col_end_index() {
    return get_parameters().polynomial_degree_plus_one_squared_ *
           std::distance(get_parameters().et_root_->begin(), cluster2_->end());
  }
  const BlockClusterTreeParameters<Scalar> &get_parameters() const {
    return *parameters_;
  }
  BlockClusterTreeParameters<Scalar> &get_parameters() { return *parameters_; }

  const TreeLeaf<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>
      &get_leaf() const {
    return *leaf_;
  }
  TreeLeaf<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> &get_leaf() {
    return *leaf_;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// private members
  //////////////////////////////////////////////////////////////////////////////
 private:
  /**
   *  \brief determines admissible clusters
   **/
  int compareCluster(const ElementTreeNode &cluster1,
                     const ElementTreeNode &cluster2) {
    int r;
    double dist = (cluster1.midpoint_ - cluster2.midpoint_).norm() -
                  cluster1.radius_ - cluster2.radius_;
    dist = dist > 0 ? dist : 0;
    double max_radius = cluster1.radius_ > cluster2.radius_ ? cluster1.radius_
                                                            : cluster2.radius_;

    /**
     * \todo replace by 2*max_radius
     */
    if (dist < 0 || max_radius >= parameters_->eta_ * dist)
      r = BlockClusterAdmissibility::Refine;
    else
      r = BlockClusterAdmissibility::LowRank;

    if (r == BlockClusterAdmissibility::Refine &&
        parameters_->max_level_ - cluster1.level_ <=
            parameters_->min_cluster_level_)
      r = BlockClusterAdmissibility::Dense;

    return r;
  }
  void updateLeafPointers() {
    // make sure we have the root of the block cluster tree calling this method
    if (cluster1_->id_ != -1 || cluster1_ != cluster2_) {
      return;
    } else {
      leaf_pointers_->clear();
      updateLeafPointers_recursion(*this);
    }
  }
  void updateLeafPointers_recursion(BlockClusterTree &child) {
    if (child.cc_ == BlockClusterAdmissibility::Refine) {
      for (auto j = 0; j < child.sons_.cols(); ++j)
        for (auto i = 0; i < child.sons_.rows(); ++i)
          updateLeafPointers_recursion(child.sons_(i, j));
    } else {
      leaf_pointers_->push_back(std::addressof(child));
    }
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// member variables
  //////////////////////////////////////////////////////////////////////////////
  std::shared_ptr<BlockClusterTreeParameters<Scalar>> parameters_;
  GenericMatrix<BlockClusterTree> sons_;
  std::shared_ptr<
      TreeLeaf<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>>
      leaf_;
  std::shared_ptr<std::vector<BlockClusterTree *>> leaf_pointers_;
  const ElementTreeNode *cluster1_;
  const ElementTreeNode *cluster2_;
  int rows_;
  int cols_;
  int cc_;
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_H2MATRIX_BLOCKCLUSTERTREE_HPP_
