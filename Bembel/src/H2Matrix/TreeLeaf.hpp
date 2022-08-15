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
//
#ifndef BEMBEL_SRC_H2MATRIX_TREELEAF_HPP_
#define BEMBEL_SRC_H2MATRIX_TREELEAF_HPP_

namespace Bembel {
template <typename Derived>
class TreeLeaf {
 public:
  //////////////////////////////////////////////////////////////////////////////
  //    Constructors
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief void constructor
   **/
  TreeLeaf(void)
      : F_(Derived(0, 0)),
        L_(Derived(0, 0)),
        R_(Derived(0, 0)),
        is_low_rank_(false) {}
  /**
   * \brief copy constructor
   **/
  TreeLeaf(const TreeLeaf &other)
      : F_(other.F_),
        L_(other.L_),
        R_(other.R_),
        is_low_rank_(other.is_low_rank_) {}
  /**
   * \brief move constructor
   **/
  TreeLeaf(TreeLeaf &&other)
      : F_(std::move(other.F_)),
        L_(std::move(other.L_)),
        R_(std::move(other.R_)),
        is_low_rank_(other.is_low_rank_) {}
  /**
   * \brief lowRank constructor
   *        whatever Eigen object is put in here will be evaluated
   */
  template <typename otherDerived>
  TreeLeaf(const Eigen::MatrixBase<otherDerived> &L,
           const Eigen::MatrixBase<otherDerived> &R)
      : F_(Derived(0, 0)), L_(L), R_(R), is_low_rank_(true) {}
  /**
   * \brief full constructor
   *        whatever Eigen object is put in here will be evaluated
   **/
  template <typename otherDerived>
  explicit TreeLeaf(const Eigen::MatrixBase<otherDerived> &F)
      : F_(F), L_(Derived(0, 0)), R_(Derived(0, 0)), is_low_rank_(false) {}
  /**
   * \brief lowRank move constructor
   **/
  TreeLeaf(Derived &&L, Derived &&R)
      : F_(Derived(0, 0)), L_(L), R_(R), is_low_rank_(true) {}
  /**
   * \brief full move constructor
   **/
  explicit TreeLeaf(Derived &&F)
      : F_(F), L_(Derived(0, 0)), R_(Derived(0, 0)), is_low_rank_(false) {}
  //////////////////////////////////////////////////////////////////////////////
  //    getter
  //////////////////////////////////////////////////////////////////////////////
  bool is_low_rank() { return is_low_rank_; }
  Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
      &get_F() {
    return F_;
  }
  Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
      &get_L() {
    return L_;
  }
  Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
      &get_R() {
    return R_;
  }
  const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
      &get_F() const {
    return F_;
  }
  const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
      &get_L() const {
    return L_;
  }
  Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
      &get_R() const {
    return R_;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    setter
  //////////////////////////////////////////////////////////////////////////////
  void set_low_rank_flag(bool flag) {
    is_low_rank_ = flag;
    return;
  }
  void set_F(const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic,
                                 Eigen::Dynamic> &F) {
    F_ = F;
  }
  void set_L(const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic,
                                 Eigen::Dynamic> &L) {
    L_ = L;
  }
  void set_R(const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic,
                                 Eigen::Dynamic> &R) {
    R_ = R;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    Operators
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief assignment operator, works for copy and move assignment
   **/
  TreeLeaf &operator=(TreeLeaf other) {
    F_.swap(other.F_);
    L_.swap(other.L_);
    R_.swap(other.R_);
    std::swap(is_low_rank_, other.is_low_rank_);
    return *this;
  }

  //////////////////////////////////////////////////////////////////////////////
  /// private members
  //////////////////////////////////////////////////////////////////////////////
 private:
  // independently how the leave is instatiated, it will always be casted to
  // an actual matrix
  Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> F_,
      L_, R_;
  bool is_low_rank_;
};
}  // namespace Bembel
#endif  // BEMBEL_SRC_H2MATRIX_TREELEAF_HPP_
