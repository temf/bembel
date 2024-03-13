// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2024 see <http://www.bembel.eu>
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
 * \ingroup AnsatzSpace
 * \brief The AnsatzSpace is the class that handles the assembly of the
 *discrete basis.
 *
 * It invokes a SuperSpace and uses the Glue and Projector class to
 *assemble a transformation matrix, which relates the SuperSpace to the desired
 *basis.
 */
template <typename Derived>
class AnsatzSpace {
 public:
  enum { Form = LinearOperatorTraits<Derived>::Form };
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Default constructor
   */
  AnsatzSpace() {}
  /**
   * \brief Copy constructor
   * \param other The object to copy from
   */
  AnsatzSpace(const AnsatzSpace &other) {
    super_space_ = other.super_space_;
    knot_repetition_ = other.knot_repetition_;
    transformation_matrix_ = other.transformation_matrix_;
  }
  /**
   * @brief Move constructor
   * @param other The object to move from
   */
  AnsatzSpace(AnsatzSpace &&other) {
    super_space_ = other.super_space_;
    knot_repetition_ = other.knot_repetition_;
    transformation_matrix_ = other.transformation_matrix_;
  }
  /**
   * \brief Copy assignment operator.
   *
   * This operator assigns the contents of another AnsatzSpace object to this
   * one.
   *
   * \param other The AnsatzSpace object to copy from.
   * \return A reference to the updated AnsatzSpace object.
   */
  AnsatzSpace &operator=(AnsatzSpace other) {
    super_space_ = other.super_space_;
    knot_repetition_ = other.knot_repetition_;
    transformation_matrix_ = other.transformation_matrix_;
    return *this;
  }
  /**
   * \brief Constructor for AnsatzSpace.
   *
   * This constructor initializes an AnsatzSpace object with the provided
   * parameters.
   *
   * \param geometry The geometry object defining the space.
   * \param refinement_level The refinement level of the space.
   * \param polynomial_degree The degree of polynomials used in the space.
   * \param knot_repetition (optional) The number of repetitions of knots in the
   * space.
   */
  AnsatzSpace(const Geometry &geometry, int refinement_level,
              int polynomial_degree, int knot_repetition = 1) {
    init_AnsatzSpace(geometry, refinement_level, polynomial_degree,
                     knot_repetition);
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_Ansatzspace
  //////////////////////////////////////////////////////////////////////////////
  /**
   * \brief Initializes the AnsatzSpace object.
   *
   * This function initializes the AnsatzSpace object with the provided
   * parameters. It sets the knot repetition, initializes the super space,
   * utilizes the Projector and Glue class to create a transformation matrix
   * which assembles conforming B-Splines from local Bernstein polynomials.
   *
   * \param geometry The geometry object.
   * \param refinement_level The refinement level of the space.
   * \param polynomial_degree The degree of polynomials used in the space.
   * \param knot_repetition The number of repetitions of knots in the space.
   */
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
  /**
   * \brief Retrieves the reference to the SuperSpace associated with this
   * AnsatzSpace.
   *
   * \return A const reference to the SuperSpace.
   */
  const SuperSpace<Derived> &get_superspace() const { return super_space_; }

  /**
   * \brief Retrieves the knot repetition value of this AnsatzSpace.
   *
   * \return The knot repetition value.
   */
  int get_knot_repetition() const { return knot_repetition_; }

  /**
   * \brief Retrieves the refinement level of this AnsatzSpace.
   *
   * \return The refinement level.
   */
  int get_refinement_level() const {
    return super_space_.get_refinement_level();
  }

  /**
   * \brief Retrieves the polynomial degree of this AnsatzSpace.
   *
   * \return The polynomial degree.
   */
  int get_polynomial_degree() const {
    return super_space_.get_polynomial_degree();
  }

  /**
   * \brief Retrieves the number of elements of the underlying ElementTree in
   * the SuperSpace of this AnsatzSpace.
   *
   * \return The number of elements.
   */
  int get_number_of_elements() const {
    return super_space_.get_number_of_elements();
  }

  /**
   * \brief Retrieves the number of patches of the underlying Geometry.
   *
   * \return The number of patches.
   */
  int get_number_of_patches() const {
    return super_space_.get_number_of_patches();
  }

  /**
   * \brief Retrieves the number of degrees of freedom of this AnsatzSpace.
   *
   * \return The number of degrees of freedom.
   */
  int get_number_of_dofs() const { return transformation_matrix_.cols(); }

  /**
   * \brief Retrieves the geometry associated with this AnsatzSpace.
   *
   * \return A const reference to the geometry.
   */
  const PatchVector &get_geometry() const {
    return super_space_.get_geometry();
  }

  /**
   * \brief Retrieves the transformation matrix associated with this
   * AnsatzSpace.
   *
   * \return A const reference to the transformation matrix.
   */
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
