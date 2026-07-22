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
#ifndef BEMBEL_SRC_SPLINE_SHAPEFUNCTIONS_HPP_
#define BEMBEL_SRC_SPLINE_SHAPEFUNCTIONS_HPP_

namespace Bembel {
namespace Basis {

/**
 *  \ingroup Spline
 *  \brief These routines implement a template recursion that allows to choose a
 *compile time instantiation of a basis-evaluation routine with a runtime p. To
 *replace the underlying basis, only these routines should be changed.
 **/
template <int P>
class PSpecificShapeFunctionHandler {
 public:
  /**
   * \brief Evaluates the Bernstein basis of degree p at point x into ar.
   */
  inline static void evalBasis(const int p, double* ar, const double x) {
    return p == P ? Bembel::Basis::EvalBernsteinBasis<P>(ar, x)
                  : PSpecificShapeFunctionHandler<P - 1>::evalBasis(p, ar, x);
  }
  /**
   * \brief Evaluates the Bernstein basis of degree p at point x into an
   * Eigen::VectorXd.
   *
   * \attention For performance critical applications use the pointer version of
   * this function to avoid slow memory allocation.
   */
  inline static Eigen::Matrix<double, Eigen::Dynamic, 1> evalBasis(
      const int p, const double x) {
    Eigen::VectorXd eval(p + 1);
    evalBasis(p, eval.data(), x);
    return eval;
  }
  /**
   * \brief Evaluates the derivative of the Bernstein basis of degree p at point
   * x into ar.
   */
  inline static void evalDerBasis(const int p, double* ar, const double x) {
    return p == P
               ? Bembel::Basis::EvalBernsteinDerBasis<P>(ar, x)
               : PSpecificShapeFunctionHandler<P - 1>::evalDerBasis(p, ar, x);
  }
  /**
   * \brief Evaluates the derivative of the Bernstein basis of degree p into an
   * Eigen::VectorXd.
   *
   * \attention For performance critical applications use the pointer version of
   * this function to avoid slow memory allocation.
   */
  inline static Eigen::Matrix<double, Eigen::Dynamic, 1> evalDerBasis(
      const int p, const double x) {
    Eigen::VectorXd eval(p + 1);
    evalDerBasis(p, eval.data(), x);
    return eval;
  }
};

template <>
class PSpecificShapeFunctionHandler<0> {
 public:
  /**
   * \brief Evaluates the Bernstein basis of degree p at point x into ar.
   */
  inline static void evalBasis(const int p, double* ar, const double x) {
    Bembel::Basis::EvalBernsteinBasis<0>(ar, x);
    return;
  }
  /**
   * \brief Evaluates the Bernstein basis of degree p at point x into an
   * Eigen::VectorXd.
   *
   * \attention For performance critical applications use the pointer version of
   * this function to avoid slow memory allocation.
   */
  inline static Eigen::VectorXd evalBasis(const int p, const double x) {
    Eigen::VectorXd eval(p + 1);
    evalBasis(p, eval.data(), x);
    return eval;
  }
  /**
   * \brief Evaluates the derivative of the Bernstein basis of degree p at point
   * x into ar.
   */
  inline static void evalDerBasis(const int p, double* ar, const double x) {
    Bembel::Basis::EvalBernsteinDerBasis<0>(ar, x);
    return;
  }
  /**
   * \brief Evaluates the derivative of the Bernstein basis of degree p into an
   * Eigen::VectorXd.
   *
   * \attention For performance critical applications use the pointer version of
   * this function to avoid slow memory allocation.
   */
  inline static Eigen::VectorXd evalDerBasis(const int p, const double x) {
    Eigen::VectorXd eval(p + 1);
    evalDerBasis(p, eval.data(), x);
    return eval;
  }
};

using ShapeFunctionHandler = PSpecificShapeFunctionHandler<Constants::MaxP>;

}  // namespace Basis
}  // namespace Bembel
#endif  // BEMBEL_SRC_SPLINE_SHAPEFUNCTIONS_HPP_
