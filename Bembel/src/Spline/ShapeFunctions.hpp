// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
//
#ifndef BEMBEL_SPLINE_SHAPEFUNCTIONS_H_
#define BEMBEL_SPLINE_SHAPEFUNCTIONS_H_

namespace Bembel {
namespace Basis {

using funptr_doubleOut_doubleptrDoubleIn = double (*)(double*, double);
using funptr_voidOut_doubleptrDoubleIn = void (*)(double*, double);

/**
 *  \ingroup Spline
 *  \brief These routines implement a template recursion that allows to choose a
 *compile time instantiation of a basis-evaluation routine with a runtime p. To
 *replace the underlying basis, only these routines should be changed.
 **/
template <int P>
class PSpecificShapeFunctionHandler {
 public:
  inline static double evalCoef(int p, double* ar, double x) {
    return p == P ? Bembel::Basis::EvalBernstein<double, P>(ar, x)
                  : PSpecificShapeFunctionHandler<P - 1>::evalCoef(p, ar, x);
  }
  inline static double evalDerCoef(int p, double* ar, double x) {
    return p == P ? Bembel::Basis::EvalBernsteinDer<double, P>(ar, x)
                  : PSpecificShapeFunctionHandler<P - 1>::evalDerCoef(p, ar, x);
  }
  inline static void evalBasis(int p, double* ar, double x) {
    return p == P ? Bembel::Basis::EvalBernsteinBasis<double, P>(ar, x)
                  : PSpecificShapeFunctionHandler<P - 1>::evalBasis(p, ar, x);
  }
  inline static void evalDerBasis(int p, double* ar, double x) {
    return p == P
               ? Bembel::Basis::EvalBernsteinDerBasis<double, P>(ar, x)
               : PSpecificShapeFunctionHandler<P - 1>::evalDerBasis(p, ar, x);
  }
  inline static constexpr funptr_doubleOut_doubleptrDoubleIn ptrEvalCoef(
      int p) {
    return p == P ? &Bembel::Basis::EvalBernstein<double, P>
                  : PSpecificShapeFunctionHandler<P - 1>::ptrEvalCoef(p);
  }
  inline static constexpr funptr_doubleOut_doubleptrDoubleIn ptrEvalDerCoef(
      int p) {
    return p == P ? &Bembel::Basis::EvalBernsteinDer<double, P>
                  : PSpecificShapeFunctionHandler<P - 1>::ptrEvalDerCoef(p);
  }
  inline static constexpr funptr_voidOut_doubleptrDoubleIn ptrEvalBasis(int p) {
    return p == P ? &Bembel::Basis::EvalBernsteinBasis<double, P>
                  : PSpecificShapeFunctionHandler<P - 1>::ptrEvalBasis(p);
  }
  inline static constexpr funptr_voidOut_doubleptrDoubleIn ptrEvalDerBasis(
      int p) {
    return p == P ? &Bembel::Basis::EvalBernsteinDerBasis<double, P>
                  : PSpecificShapeFunctionHandler<P - 1>::ptrEvalDerBasis(p);
  }
  inline static constexpr bool checkP(int p) {
    static_assert(P > 0, "Polynomial degree must be larger than zero");
    return p <= Constants::MaxP;
  }
};

template <>
class PSpecificShapeFunctionHandler<0> {
 public:
  inline static double evalCoef(int p, double* ar, double x) {
    return Bembel::Basis::EvalBernstein<double, 0>(ar, x);
  }
  inline static double evalDerCoef(int p, double* ar, double x) {
    return Bembel::Basis::EvalBernsteinDer<double, 0>(ar, x);
  }
  inline static void evalBasis(int p, double* ar, double x) {
    return Bembel::Basis::EvalBernsteinBasis<double, 0>(ar, x);
  }
  inline static void evalDerBasis(int p, double* ar, double x) {
    return Bembel::Basis::EvalBernsteinDerBasis<double, 0>(ar, x);
  }
  inline static constexpr funptr_doubleOut_doubleptrDoubleIn ptrEvalCoef(
      int p) {
    return &Bembel::Basis::EvalBernstein<double, 0>;
  }
  inline static constexpr funptr_doubleOut_doubleptrDoubleIn ptrEvalDerCoef(
      int p) {
    return &Bembel::Basis::EvalBernsteinDer<double, 0>;
  }
  inline static constexpr funptr_voidOut_doubleptrDoubleIn ptrEvalBasis(int p) {
    return &Bembel::Basis::EvalBernsteinBasis<double, 0>;
  }
  inline static constexpr funptr_voidOut_doubleptrDoubleIn ptrEvalDerBasis(
      int p) {
    return &Bembel::Basis::EvalBernsteinDerBasis<double, 0>;
  }
  inline static constexpr bool checkP(int p) { return Constants::MaxP >= 0; }
};

using ShapeFunctionHandler = PSpecificShapeFunctionHandler<Constants::MaxP>;

}  // namespace Basis
}  // namespace Bembel
#endif
