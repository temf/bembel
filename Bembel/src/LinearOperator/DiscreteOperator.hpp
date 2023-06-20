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
#ifndef BEMBEL_SRC_LINEAROPERATOR_DISCRETEOPERATOR_HPP_
#define BEMBEL_SRC_LINEAROPERATOR_DISCRETEOPERATOR_HPP_

namespace Bembel {
/**
 *  \ingroup LinearOperator
 *  \brief Helper struct that is used in order to partially specialise the
 *         compute routine of DiscreteOperator for different types of
 *         output formats
 */
template <typename MatrixFormat, typename Derived>
struct DiscreteOperatorComputer {};
/**
 *  \brief Helper struct that is used in order to partially specialise the
 *         compute routine of DiscreteOperator for the Eigen::MatrixXd format
 */
template <typename Derived, typename Scalar>
struct DiscreteOperatorComputer<
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>, Derived> {
  DiscreteOperatorComputer() {}
  static void compute(
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> *disc_op,
      const Derived &lin_op, const AnsatzSpace<Derived> &ansatz_space) {
    GaussSquare<Constants::maximum_quadrature_degree> GS;
    const SuperSpace<Derived> &super_space = ansatz_space.get_superspace();
    const ElementTree &element_tree = super_space.get_mesh().get_element_tree();
    auto number_of_elements = element_tree.get_number_of_elements();
    const auto vector_dimension =
        getFunctionSpaceVectorDimension<LinearOperatorTraits<Derived>::Form>();
    auto polynomial_degree = super_space.get_polynomial_degree();
    auto polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);
    auto ffield_deg = lin_op.get_FarfieldQuadratureDegree(polynomial_degree);
    auto ffield_qnodes =
        DuffyTrick::computeFfieldQnodes(super_space, GS[ffield_deg]);
    disc_op->resize(vector_dimension * polynomial_degree_plus_one_squared *
                        number_of_elements,
                    vector_dimension * polynomial_degree_plus_one_squared *
                        number_of_elements);
#pragma omp parallel
    {
#pragma omp single
      {
        for (auto element1 = element_tree.cpbegin();
             element1 != element_tree.cpend(); ++element1)
          for (auto element2 = element_tree.cpbegin();
               element2 != element_tree.cpend(); ++element2) {
#pragma omp task firstprivate(element1, element2)
            {
              Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> intval(
                  vector_dimension * polynomial_degree_plus_one_squared,
                  vector_dimension * polynomial_degree_plus_one_squared);
              DuffyTrick::evaluateBilinearForm(
                  lin_op, super_space, *element1, *element2, GS,
                  ffield_qnodes[element1->id_], ffield_qnodes[element2->id_],
                  &intval);
              for (auto i = 0; i < vector_dimension; ++i)
                for (auto j = 0; j < vector_dimension; ++j)
                  disc_op->block(polynomial_degree_plus_one_squared *
                                     (i * number_of_elements + element1->id_),
                                 polynomial_degree_plus_one_squared *
                                     (j * number_of_elements + element2->id_),
                                 polynomial_degree_plus_one_squared,
                                 polynomial_degree_plus_one_squared) =
                      intval.block(i * polynomial_degree_plus_one_squared,
                                   j * polynomial_degree_plus_one_squared,
                                   polynomial_degree_plus_one_squared,
                                   polynomial_degree_plus_one_squared);
            }
          }
      }
    }
    auto projector = ansatz_space.get_transformation_matrix();
    disc_op[0] = projector.transpose() * (disc_op[0] * projector);
    return;
  }
};
/**
 *  \brief Helper struct that is used in order to partially specialise the
 *         compute routine of DiscreteOperator for the Eigen::H2Matrix format
 */
template <typename Derived>
struct DiscreteOperatorComputer<
    Eigen::H2Matrix<typename LinearOperatorTraits<Derived>::Scalar>, Derived> {
  DiscreteOperatorComputer() {}
  static void compute(
      Eigen::H2Matrix<typename LinearOperatorTraits<Derived>::Scalar> *disc_op,
      const Derived &lin_op, const AnsatzSpace<Derived> &ansatz_space) {
    disc_op->init_H2Matrix(lin_op, ansatz_space);
    return;
  }
};
/**
 *  \brief DiscreteOperator
 *  \todo  Add a documentation
 */
template <typename MatrixFormat, typename Derived>
class DiscreteOperator {
 public:
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  DiscreteOperator() {}
  explicit DiscreteOperator(const AnsatzSpace<Derived> &ansatz_space) {
    init_DiscreteOperator(ansatz_space);
  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_DiscreteOperator
  //////////////////////////////////////////////////////////////////////////////
  void init_DiscreteOperator(const AnsatzSpace<Derived> &ansatz_space) {
    ansatz_space_ = ansatz_space;
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    compute
  //////////////////////////////////////////////////////////////////////////////
  inline void compute() {
    DiscreteOperatorComputer<MatrixFormat, Derived>::compute(&disc_op_, lin_op_,
                                                             ansatz_space_);
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    getter
  //////////////////////////////////////////////////////////////////////////////
  const MatrixFormat &get_discrete_operator() const { return disc_op_; }
  const Derived &get_linear_operator() const { return lin_op_; }
  Derived &get_linear_operator() { return lin_op_; }
  //////////////////////////////////////////////////////////////////////////////
  //    private member variables
  //////////////////////////////////////////////////////////////////////////////
 private:
  Derived lin_op_;
  MatrixFormat disc_op_;
  AnsatzSpace<Derived> ansatz_space_;
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_LINEAROPERATOR_DISCRETEOPERATOR_HPP_
