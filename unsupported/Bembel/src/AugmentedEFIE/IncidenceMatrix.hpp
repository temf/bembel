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
#ifndef BEMBEL_SRC_AUGMENTEDEFIE_INCIDENCEMATRIX_HPP_
#define BEMBEL_SRC_AUGMENTEDEFIE_INCIDENCEMATRIX_HPP_

namespace Bembel {
/**
 *  \ingroup AugmentedEFIE
 */
template <typename Derived_vec, typename Derived>
class IncidenceMatrix {
  typedef Eigen::SparseMatrix<typename LinearOperatorTraits<Derived>::Scalar>
      MatrixType;
  typedef Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                        Eigen::Dynamic, 1>
      VectorType;

 public:
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  IncidenceMatrix() {}
  IncidenceMatrix(const AnsatzSpace<Derived_vec> &ansatz_space_vec,
                  const AnsatzSpace<Derived> &ansatz_space) {
    static_assert(
        std::is_same<
            typename LinearOperatorTraits<Derived>::Scalar,
            typename LinearOperatorTraits<Derived_vec>::Scalar>::value &&
        "Scalar Types need to be equal");
    init_IncidenceMatrix(ansatz_space_vec, ansatz_space);
  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_IncidenceMatrix
  //////////////////////////////////////////////////////////////////////////////
  void init_IncidenceMatrix(const AnsatzSpace<Derived_vec> &ansatz_space_vec,
                            const AnsatzSpace<Derived> &ansatz_space) {
    ansatz_space_vec_ = ansatz_space_vec;
    ansatz_space_ = ansatz_space;
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    compute
  //////////////////////////////////////////////////////////////////////////////
  void compute() {
    const SuperSpace<Derived_vec> &super_space_vec =
        ansatz_space_vec_.get_superspace();
    const SuperSpace<Derived> &super_space = ansatz_space_.get_superspace();
    const auto &element_tree = super_space.get_mesh().get_element_tree();

    // Quadrature
    GaussSquare<20> GS;
    Cubature Q = GS[19];
    const auto &number_of_elements = element_tree.get_number_of_elements();
    const auto polynomial_degree_plus_one_squared_vec =
        super_space_vec.get_polynomial_degree_plus_one_squared();
    const auto polynomial_degree_plus_one_squared =
        super_space.get_polynomial_degree_plus_one_squared();

    // Triplets
    typedef Eigen::Triplet<typename LinearOperatorTraits<Derived>::Scalar> T;
    std::vector<T> tripletList;
    tripletList.reserve(2 * polynomial_degree_plus_one_squared_vec *
                        polynomial_degree_plus_one_squared *
                        number_of_elements);

    // Iterate over elements
    for (auto element = element_tree.cpbegin(); element != element_tree.cpend();
         ++element) {
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, Eigen::Dynamic>
          intval(2 * polynomial_degree_plus_one_squared_vec,
                 polynomial_degree_plus_one_squared);
      intval.setZero();
      // Iterate over quadrature points
      for (auto i = 0; i < Q.w_.size(); ++i) {
        SurfacePoint qp;
        super_space.map2surface(*element, Q.xi_.col(i), Q.w_(i), &qp);
        evaluateIntegrand(super_space_vec, super_space, qp, &intval);
      }
      // Transform local element matrices to triplets
      for (auto jj = 0; jj < polynomial_degree_plus_one_squared; ++jj)
        for (auto i = 0; i < 2; ++i)
          for (auto ii = 0; ii < polynomial_degree_plus_one_squared_vec; ++ii)
            tripletList.push_back(
                T(polynomial_degree_plus_one_squared_vec *
                          (i * number_of_elements + element->id_) +
                      ii,
                  polynomial_degree_plus_one_squared * element->id_ + jj,
                  intval(i * polynomial_degree_plus_one_squared_vec + ii, jj)));
    }
    incidence_matrix_.resize(
        2 * polynomial_degree_plus_one_squared_vec * number_of_elements,
        polynomial_degree_plus_one_squared * number_of_elements);
    incidence_matrix_.setFromTriplets(tripletList.begin(), tripletList.end());
    auto projector_vec = ansatz_space_vec_.get_transformation_matrix();
    auto projector = ansatz_space_.get_transformation_matrix();
    incidence_matrix_ =
        projector_vec.transpose() * (incidence_matrix_ * projector);
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    getter
  //////////////////////////////////////////////////////////////////////////////
  const MatrixType &get_incidence_matrix() const { return incidence_matrix_; }
  //////////////////////////////////////////////////////////////////////////////
  //    private member variables
  //////////////////////////////////////////////////////////////////////////////
 private:
  template <class T_vec, class T>
  void evaluateIntegrand(
      const T_vec &super_space_vec, const T &super_space, const SurfacePoint &p,
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, Eigen::Dynamic> *intval) {
    const auto polynomial_degree_plus_one_squared_vec =
        super_space_vec.get_polynomial_degree_plus_one_squared();

    // const auto polynomial_degree_plus_one_squared =
    //     super_space.get_polynomial_degree_plus_one_squared();
    VectorType intval_scal =
        super_space.basis(p.segment<2>(0)).array().real();
    VectorType intval_dx = super_space_vec.basisDx(p.segment<2>(0));
    VectorType intval_dy = super_space_vec.basisDy(p.segment<2>(0));
    VectorType intval_div(2 * polynomial_degree_plus_one_squared_vec);
    intval_div << intval_dx, intval_dy;

    // get basic information
    const unsigned int elements_per_direction =
        (1 << super_space.get_refinement_level());
    const double h = 1. / ((double)(elements_per_direction));

    // compute surface measures from tangential derivatives
    double x_kappa = p.segment<3>(6).cross(p.segment<3>(9)).norm();
    double w = p(2) / h;
    intval[0] += w * intval_div * intval_scal.transpose();
  }
  MatrixType incidence_matrix_;
  AnsatzSpace<Derived_vec> ansatz_space_vec_;
  AnsatzSpace<Derived> ansatz_space_;
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_AUGMENTEDEFIE_INCIDENCEMATRIX_HPP_
