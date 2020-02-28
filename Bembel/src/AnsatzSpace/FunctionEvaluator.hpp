// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_ANSATZSPACE_FUNCTIONEVALUATOR_H_
#define BEMBEL_ANSATZSPACE_FUNCTIONEVALUATOR_H_

namespace Bembel {
/**
 *  \ingroup AnsatzSpace
 *  \brief The FunctionEvaluator provides means to evaluate coefficient vectors
 * as functions on the geometry.
 */
template <typename Derived>
class FunctionEvaluator {
 public:
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  FunctionEvaluator() {}
  FunctionEvaluator(const FunctionEvaluator &other) = default;
  FunctionEvaluator(FunctionEvaluator &&other) = default;
  FunctionEvaluator &operator=(FunctionEvaluator other) {
    ansatz_space_ = other.ansatz_space_;
    fun_ = other.fun_;
    polynomial_degree_plus_one_squared_ =
        other.polynomial_degree_plus_one_squared_;
    return *this;
  }
  FunctionEvaluator(const AnsatzSpace<Derived> &ansatz_space) {
    init_FunctionEvaluator(ansatz_space);
    return;
  }
  FunctionEvaluator(
      const AnsatzSpace<Derived> &ansatz_space,
      const Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                          Eigen::Dynamic, 1> &fun) {
    init_FunctionEvaluator(ansatz_space, fun);
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_Ansatzspace
  //////////////////////////////////////////////////////////////////////////////
  void init_FunctionEvaluator(const AnsatzSpace<Derived> &ansatz_space) {
    ansatz_space_ = ansatz_space;
    auto polynomial_degree = ansatz_space_.get_polynomial_degree();
    polynomial_degree_plus_one_squared_ =
        (polynomial_degree + 1) * (polynomial_degree + 1);
    reordering_vector_ = ansatz_space_.get_superspace()
                             .get_mesh()
                             .get_element_tree()
                             .computeReorderingVector();
    return;
  }
  void init_FunctionEvaluator(
      const AnsatzSpace<Derived> &ansatz_space,
      const Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                          Eigen::Dynamic, 1> &fun) {
    ansatz_space_ = ansatz_space;
    set_function(fun);
    auto polynomial_degree = ansatz_space_.get_polynomial_degree();
    polynomial_degree_plus_one_squared_ =
        (polynomial_degree + 1) * (polynomial_degree + 1);
    reordering_vector_ = ansatz_space_.get_superspace()
                             .get_mesh()
                             .get_element_tree()
                             .computeReorderingVector();
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    evaluators
  //////////////////////////////////////////////////////////////////////////////
  Eigen::Matrix<
      typename LinearOperatorTraits<Derived>::Scalar,
      getFunctionSpaceOutputDimension<LinearOperatorTraits<Derived>::Form>(), 1>
  evaluateOnPatch(int patch, const Eigen::Vector2d &ref_point) const {
    const int elements_per_direction =
        (1 << ansatz_space_.get_refinement_level());
    const int elements_per_patch =
        elements_per_direction * elements_per_direction;
    const double h = 1. / ((double)(elements_per_direction));
    const int x_idx_ = std::floor(ref_point(0) / h);
    const int y_idx_ = std::floor(ref_point(1) / h);
    const int x_idx = std::min(std::max(x_idx_, 0), elements_per_direction - 1);
    const int y_idx = std::min(std::max(y_idx_, 0), elements_per_direction - 1);
    const int tp_idx =
        x_idx + elements_per_direction * y_idx + patch * elements_per_patch;
    const int et_idx = reordering_vector_[tp_idx];
    const ElementTreeNode &element = *(
        ansatz_space_.get_superspace().get_mesh().get_element_tree().cpbegin() +
        et_idx);

    SurfacePoint sp;
    ansatz_space_.get_superspace().get_geometry()[patch].updateSurfacePoint(
        &sp, ref_point, 1, element.mapToReferenceElement(ref_point));
    return evaluate(element, sp);
  }
  Eigen::Matrix<
      typename LinearOperatorTraits<Derived>::Scalar,
      getFunctionSpaceOutputDimension<LinearOperatorTraits<Derived>::Form>(), 1>
  evaluate(const ElementTreeNode &element, const SurfacePoint &p) const {
    return eval_.eval(
        ansatz_space_.get_superspace(), polynomial_degree_plus_one_squared_,
        element, p,
        fun_.block(polynomial_degree_plus_one_squared_ * element.id_, 0,
                   polynomial_degree_plus_one_squared_,
                   getFunctionSpaceVectorDimension<
                       LinearOperatorTraits<Derived>::Form>()));
  }
  typename LinearOperatorTraits<Derived>::Scalar evaluateDiv(
      const ElementTreeNode &element, const SurfacePoint &p) const {
    return eval_.evalDiv(
        ansatz_space_.get_superspace(), polynomial_degree_plus_one_squared_,
        element, p,
        fun_.block(polynomial_degree_plus_one_squared_ * element.id_, 0,
                   polynomial_degree_plus_one_squared_,
                   getFunctionSpaceVectorDimension<
                       LinearOperatorTraits<Derived>::Form>()));
  }
  //////////////////////////////////////////////////////////////////////////////
  //    setters
  //////////////////////////////////////////////////////////////////////////////
  void set_function(
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, 1>
          fun) {
    const auto vec_dim =
        getFunctionSpaceVectorDimension<LinearOperatorTraits<Derived>::Form>();
    auto longfun = (ansatz_space_.get_transformation_matrix() * fun).eval();
    fun_ =
        Eigen::Map<Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                                 Eigen::Dynamic, vec_dim>>(
            longfun.data(), longfun.rows() / vec_dim, vec_dim);
  }
  //////////////////////////////////////////////////////////////////////////////
  //    private member variables
  //////////////////////////////////////////////////////////////////////////////
 private:
  std::vector<int> reordering_vector_;
  AnsatzSpace<Derived> ansatz_space_;
  Eigen::Matrix<
      typename LinearOperatorTraits<Derived>::Scalar, Eigen::Dynamic,
      getFunctionSpaceVectorDimension<LinearOperatorTraits<Derived>::Form>()>
      fun_;
  int polynomial_degree_plus_one_squared_;
  FunctionEvaluatorEval<typename LinearOperatorTraits<Derived>::Scalar,
                        LinearOperatorTraits<Derived>::Form, Derived>
      eval_;
};

}  // namespace Bembel
#endif
