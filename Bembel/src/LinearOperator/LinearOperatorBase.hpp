// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_LINEAROPERATOR_LINEAROPERATORBASE_H_
#define BEMBEL_LINEAROPERATOR_LINEAROPERATORBASE_H_

namespace Bembel {
/**
 *    \ingroup LinearOperator
 *    \brief linear operator base class. this serves as a common interface for
 *           existing linear operators
 **/
template <typename Derived>
struct LinearOperatorBase {
  // Constructors
  LinearOperatorBase(){};
  // the user has to provide the implementation of this function, which
  // is able to evaluate the integrand of the Galerkin formulation in a
  // pair of quadrature points represented as a
  // Surface point [xi; h * w; Chi(xi); dsChi(xi); dtChi(xi)]
  template <class T>
  void evaluateIntegrand(
      const T &super_space, const SurfacePoint &p1, const SurfacePoint &p2,
      Eigen::Matrix<typename LinearOperatorTraits<Derived>::Scalar,
                    Eigen::Dynamic, Eigen::Dynamic> *intval) const {
    static_cast<const Derived *>(this)->evaluateIntegrand_impl(super_space, p1,
                                                               p2, intval);
    return;
  }
  // the user has to provide the implementation of this function, which
  // is able to evaluate the interpolation values of the Galerkin formulation on
  // the reference domain in a pair of interpolation points represented as a
  // Surface point [xi; 1.; Chi(xi); dsChi(xi); dtChi(xi)]
  Eigen::Matrix<
      typename LinearOperatorTraits<Derived>::Scalar,
      getFunctionSpaceVectorDimension<LinearOperatorTraits<Derived>::Form>() *
          LinearOperatorTraits<Derived>::NumberOfFMMComponents,
      getFunctionSpaceVectorDimension<LinearOperatorTraits<Derived>::Form>() *
          LinearOperatorTraits<Derived>::NumberOfFMMComponents>
  evaluateFMMInterpolation(const SurfacePoint &p1,
                           const SurfacePoint &p2) const {
    return static_cast<const Derived *>(this)->evaluateFMMInterpolation_impl(
        p1, p2);
  }
  // return the required quadrature degree for the far-field
  int get_FarfieldQuadratureDegree(int ansatz_degree) const {
    return ansatz_degree - LinearOperatorTraits<Derived>::OperatorOrder + 1;
  }
  /**
   * \brief Compute quadrature degree for numerical integretation close to the
   *        singularity based on distance, refinement level, degree of ansatz
   *        functions and operator_order.
   *
   * See also WAVELET GALERKIN SCHEMES FOR BOUNDARY INTEGRAL EQUATIONS -
   * IMPLEMENTATION AND QUADRATURE by Helmut Harbrecht and Reinhold Schneider
   * for more details.
   **/
  // return the required quadrature degree for the near-field
  int getNearfieldQuadratureDegree(int ansatz_degree, double distance,
                                   int level) const {
    // if farfield quadrature is computed: take log of distance
    // otherwise compute quadrature degree for Duffy trick
    double distance_log =
        (distance * (1 << level) < 1) ? -(level * log(2.)) : log(distance);

    // alpha/2 is the convergence rate of the corresponding Galerkin solution,
    // this ensures the correct convergence rate of the potential
    int alpha =
        2 - LinearOperatorTraits<Derived>::OperatorOrder + 2 * ansatz_degree;

    // compute numerator and denominator of quadrature degree
    double numerator =
        (alpha + ansatz_degree) * level * log(2.) -
        (2 - ansatz_degree + LinearOperatorTraits<Derived>::OperatorOrder) *
            distance_log;
    double denominator = (level + 2) * log(2.) + distance_log;

    return 0.5 * numerator / denominator;
  }
  // pointer to the derived object
  Derived &derived() { return *static_cast<Derived *>(this); }
  // const pointer to the derived object
  const Derived &derived() const { return *static_cast<const Derived *>(this); }
};
}  // namespace Bembel
#endif
