// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_POTENTIAL_DISCRETEPOTENTIAL_H_
#define BEMBEL_POTENTIAL_DISCRETEPOTENTIAL_H_

namespace Bembel {
/**
 *  \ingroup Potential
 *  \brief DiscretePotential
 *  \todo  Add a documentation
 */
template <typename Derived, typename LinOp>
class DiscretePotential {
 public:
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  DiscretePotential() {}
  DiscretePotential(const AnsatzSpace<LinOp> &ansatz_space) {
    init_DiscretePotential(ansatz_space);
  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_DiscretePotential
  //////////////////////////////////////////////////////////////////////////////
  void init_DiscretePotential(const AnsatzSpace<LinOp> &ansatz_space) {
    ansatz_space_ = ansatz_space;
    fun_ev_ = FunctionEvaluator<LinOp>(ansatz_space_);
    /**
     * \todo obtain this from ansatz space
     */
    deg_ = ansatz_space_.get_polynomial_degree() + 1;
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    compute
  //////////////////////////////////////////////////////////////////////////////
  static_assert(
      getFunctionSpaceOutputDimension<LinearOperatorTraits<LinOp>::Form>() ==
          PotentialTraits<Derived>::OutputSpaceDimension,
      "Dimension mismatch in potential evaluation");
  Eigen::Matrix<typename PotentialReturnScalar<
                    typename LinearOperatorTraits<LinOp>::Scalar,
                    typename PotentialTraits<Derived>::Scalar>::Scalar,
                Eigen::Dynamic, PotentialTraits<Derived>::OutputSpaceDimension>
  evaluate(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points) {
    auto FunctionSpaceVectorDimension =
        getFunctionSpaceVectorDimension<LinearOperatorTraits<LinOp>::Form>();
    auto OutputDimension = PotentialTraits<Derived>::OutputSpaceDimension;

    GaussSquare<Constants::maximum_quadrature_degree> GS;
    auto Q = GS[deg_];

    auto super_space = ansatz_space_.get_superspace();

    auto element_tree = super_space.get_mesh().get_element_tree();
    auto number_of_elements = element_tree.get_number_of_elements();

    auto polynomial_degree = super_space.get_polynomial_degree();
    auto polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);

    Eigen::Matrix<typename PotentialReturnScalar<
                      typename LinearOperatorTraits<LinOp>::Scalar,
                      typename PotentialTraits<Derived>::Scalar>::Scalar,
                  Eigen::Dynamic,
                  PotentialTraits<Derived>::OutputSpaceDimension>
        potential;
    potential.resize(points.rows(),
                     PotentialTraits<Derived>::OutputSpaceDimension);
    potential.setZero();

#pragma omp parallel
    {
      Eigen::Matrix<typename PotentialReturnScalar<
                        typename LinearOperatorTraits<LinOp>::Scalar,
                        typename PotentialTraits<Derived>::Scalar>::Scalar,
                    Eigen::Dynamic,
                    PotentialTraits<Derived>::OutputSpaceDimension>
          my_potential;
      my_potential.resize(points.rows(),
                          PotentialTraits<Derived>::OutputSpaceDimension);
      my_potential.setZero();
      for (auto element = element_tree.cpbegin();
           element != element_tree.cpend(); ++element) {
#pragma omp single nowait
        {
          SurfacePoint qp;
          for (auto j = 0; j < Q.w_.size(); ++j) {
            super_space.map2surface(
                *element, Q.xi_.col(j),
                element->get_h() * element->get_h() * Q.w_(j), &qp);
            for (auto i = 0; i < points.rows(); ++i) {
              my_potential.row(i) += pot_.evaluateIntegrand_impl(
                  fun_ev_, *element, points.row(i), qp);
            }
          }
        }
      }
#pragma omp critical
      potential += my_potential;
    }
    return potential;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    setter
  //////////////////////////////////////////////////////////////////////////////
  void set_cauchy_data(
      const Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar,
                          Eigen::Dynamic, 1> &cauchy_data) {
    fun_ev_.set_function(cauchy_data);
  }
  void set_degree(const int &deg) { deg_ = deg; }
  //////////////////////////////////////////////////////////////////////////////
  //    getter
  //////////////////////////////////////////////////////////////////////////////
  Derived &get_potential() { return pot_; }
  //////////////////////////////////////////////////////////////////////////////
  //    private member variables
  //////////////////////////////////////////////////////////////////////////////
 private:
  int deg_;
  Derived pot_;
  AnsatzSpace<LinOp> ansatz_space_;
  FunctionEvaluator<LinOp> fun_ev_;
};  // namespace Bembel

}  // namespace Bembel
#endif
