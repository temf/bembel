// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_LINEARFORM_DISCRETELINEARFORM_H_
#define BEMBEL_LINEARFORM_DISCRETELINEARFORM_H_

namespace Bembel {
/**
 *  \ingroup LinearForm
 *  \brief This class, given a LinearForm with defined traits, takes care of the
 * assembly of the right hand side.
 */
template <typename Derived, typename LinOp>
class DiscreteLinearForm {
 public:
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  DiscreteLinearForm() {}
  DiscreteLinearForm(const AnsatzSpace<LinOp> &ansatz_space) {
    init_DiscreteLinearForm(ansatz_space);
  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_DiscreteLinearForm
  //////////////////////////////////////////////////////////////////////////////
  void init_DiscreteLinearForm(const AnsatzSpace<LinOp> &ansatz_space) {
    ansatz_space_ = ansatz_space;
    deg_ = ansatz_space_.get_polynomial_degree() + 1;
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    compute
  //////////////////////////////////////////////////////////////////////////////
  /**
   *    \todo Add inline commentary
   **/
  void compute() {
    GaussSquare<Constants::maximum_quadrature_degree> GS;
    auto Q = GS[deg_];
    SurfacePoint qp;
    auto super_space = ansatz_space_.get_superspace();
    auto element_tree = super_space.get_mesh().get_element_tree();
    auto number_of_elements = element_tree.get_number_of_elements();
    auto polynomial_degree = super_space.get_polynomial_degree();
    auto polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);
    const auto function_space_dimension =
        getFunctionSpaceVectorDimension<LinearOperatorTraits<LinOp>::Form>();
    Eigen::Matrix<typename LinearFormTraits<Derived>::Scalar, Eigen::Dynamic,
                  function_space_dimension>
        intval(polynomial_degree_plus_one_squared, function_space_dimension);
    Eigen::Matrix<typename LinearFormTraits<Derived>::Scalar, Eigen::Dynamic,
                  function_space_dimension>
        disc_lf_matrix(polynomial_degree_plus_one_squared * number_of_elements,
                       function_space_dimension);
    disc_lf_matrix.setZero();
    for (auto element = element_tree.cpbegin(); element != element_tree.cpend();
         ++element) {
      intval.setZero();
      for (auto i = 0; i < Q.w_.size(); ++i) {
        super_space.map2surface(*element, Q.xi_.col(i),
                                element->get_h() * Q.w_(i), &qp);
        lf_.evaluateIntegrand_impl(super_space, qp, &intval);
      }
      disc_lf_matrix.block(polynomial_degree_plus_one_squared * element->id_, 0,
                           polynomial_degree_plus_one_squared,
                           function_space_dimension) = intval;
    }
    disc_lf_ =
        ansatz_space_.get_transformation_matrix().transpose() *
        Eigen::Map<Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar,
                                 Eigen::Dynamic, 1>>(
            disc_lf_matrix.data(),
            disc_lf_matrix.rows() * disc_lf_matrix.cols());
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    setter
  //////////////////////////////////////////////////////////////////////////////
  void set_degree(const int &deg) { deg_ = deg; }
  //////////////////////////////////////////////////////////////////////////////
  //    getter
  //////////////////////////////////////////////////////////////////////////////
  Derived &get_linear_form() { return lf_; }
  const Eigen::Matrix<typename LinearFormTraits<Derived>::Scalar,
                      Eigen::Dynamic, 1>
      &get_discrete_linear_form() const {
    return disc_lf_;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    private member variables
  //////////////////////////////////////////////////////////////////////////////
 private:
  int deg_;
  Derived lf_;
  Eigen::Matrix<typename LinearFormTraits<Derived>::Scalar, Eigen::Dynamic, 1>
      disc_lf_;
  AnsatzSpace<LinOp> ansatz_space_;
};

}  // namespace Bembel
#endif
