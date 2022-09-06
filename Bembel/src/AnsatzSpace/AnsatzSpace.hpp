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
#ifndef BEMBEL_SRC_ANSATZSPACE_ANSATZSPACE_HPP_
#define BEMBEL_SRC_ANSATZSPACE_ANSATZSPACE_HPP_

namespace Bembel {
/**
 *  \ingroup AnsatzSpace
 *  \brief The AnsatzSpace is the class that handles the assembly of the
 *discrete basis.
 *
 *	It invokes a superspace and uses the Glue and Projektor class to
 *assemble a transformation matrix, which relates the superspace to the desired
 *basis.
 */
template <typename Derived>
class AnsatzSpace {
 public:
  enum { Form = LinearOperatorTraits<Derived>::Form };
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  AnsatzSpace() {}
  AnsatzSpace(const AnsatzSpace &other) {
    super_space_ = other.super_space_;
    knot_repetition_ = other.knot_repetition_;
    transformation_matrix_ = other.transformation_matrix_;
  }
  AnsatzSpace(AnsatzSpace &&other) {
    super_space_ = other.super_space_;
    knot_repetition_ = other.knot_repetition_;
    transformation_matrix_ = other.transformation_matrix_;
  }
  AnsatzSpace &operator=(AnsatzSpace other) {
    super_space_ = other.super_space_;
    knot_repetition_ = other.knot_repetition_;
    transformation_matrix_ = other.transformation_matrix_;
    return *this;
  }

  AnsatzSpace(const Geometry &geometry, int refinement_level,
              int polynomial_degree, int knot_repetition = 1) {
    init_AnsatzSpace(geometry, refinement_level, polynomial_degree,
                     knot_repetition);
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_Ansatzspace
  //////////////////////////////////////////////////////////////////////////////
  void init_AnsatzSpace(const Geometry &geometry, int refinement_level,
                        int polynomial_degree, int knot_repetition) {
    knot_repetition_ = knot_repetition;
    super_space_.init_SuperSpace(geometry, refinement_level, polynomial_degree);
    Projector<Derived> proj(super_space_, knot_repetition_);
    Glue<Derived> glue(super_space_, proj);
    transformation_matrix_ =
        proj.get_projection_matrix() * glue.get_glue_matrix();
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    getters
  //////////////////////////////////////////////////////////////////////////////
  const SuperSpace<Derived> &get_superspace() const { return super_space_; }
  int get_knot_repetition() const { return knot_repetition_; }
  int get_refinement_level() const {
    return super_space_.get_refinement_level();
  }
  int get_polynomial_degree() const {
    return super_space_.get_polynomial_degree();
  }
  int get_number_of_elements() const {
    return super_space_.get_number_of_elements();
  }
  int get_number_of_patches() const {
    return super_space_.get_number_of_patches();
  }
  int get_number_of_dofs() const { return transformation_matrix_.cols(); }
  const PatchVector &get_geometry() const {
    return super_space_.get_geometry();
  }
  const Eigen::SparseMatrix<double> &get_transformation_matrix() const {
    return transformation_matrix_;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    private member variables
  //////////////////////////////////////////////////////////////////////////////
 private:
  Eigen::SparseMatrix<double> transformation_matrix_;
  SuperSpace<Derived> super_space_;
  int knot_repetition_;
};
}  // namespace Bembel
#endif  // BEMBEL_SRC_ANSATZSPACE_ANSATZSPACE_HPP_
