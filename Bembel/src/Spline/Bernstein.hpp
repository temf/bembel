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

#ifndef BEMBEL_SRC_SPLINE_BERNSTEIN_HPP_
#define BEMBEL_SRC_SPLINE_BERNSTEIN_HPP_

namespace Bembel {
namespace Basis {
/**
 *  \ingroup Spline
 *  \brief Template recursion to produce Bernstein polynomials. This is only
 *         limited by the binomial coefficient, see Pascal.hpp
 */
template <int N>
inline constexpr double BernsteinX(double evaluation_point) noexcept {
#ifdef _spline_debug_flag_
  assert((evaluation_point > -.0000001) && (evaluation_point < 1.0000001) &&
         ("Function only valid for 0 <= x <= 1!"));
#endif
  return evaluation_point * BernsteinX<N - 1>(evaluation_point);
}

template <>
inline constexpr double BernsteinX<1>(double evaluation_point) noexcept {
  return evaluation_point;
}
template <>
inline constexpr double BernsteinX<0>(double evaluation_point) noexcept {
  return 1.;
}
template <>
inline constexpr double BernsteinX<-1>(double evaluation_point) noexcept {
  return 0.;
}

template <int N, int P>
inline constexpr double Bernstein(double evaluation_point) noexcept {
  return Binomial<N, P>::value * BernsteinX<N>(evaluation_point) *
         BernsteinX<P - N>(1. - evaluation_point);
}
////////////////////////////////////////////////////////////////////////////////
/// Hidden Classes
////////////////////////////////////////////////////////////////////////////////
template <int N, int P>
class HiddenBernsteinClass {
 public:
  // think twice when replacing double* by something else. All other attempts
  // were slower
  static inline void EvalBasis(double *in,
                               const double evaluation_point) noexcept {
    in[N] = Bernstein<N, P>(evaluation_point);
    HiddenBernsteinClass<N - 1, P>::EvalBasis(in, evaluation_point);
    return;
  }
  // think twice when replacing double* by something else. All other attempts
  // were slower
  static inline void EvalDerBasis(double *in,
                                  const double evaluation_point) noexcept {
    in[N] = (P + 1) * (Bernstein<N - 1, P>(evaluation_point) -
                       Bernstein<N, P>(evaluation_point));
    HiddenBernsteinClass<N - 1, P>::EvalDerBasis(in, evaluation_point);
    return;
  }
};

template <int P>
class HiddenBernsteinClass<0, P> {
 public:
  static inline void EvalBasis(double *in,
                               const double evaluation_point) noexcept {
    in[0] = Bernstein<0, P>(evaluation_point);
    return;
  }
  static inline void EvalDerBasis(double *in,
                                  const double evaluation_point) noexcept {
    in[0] = (-P - 1) * Bernstein<0, P>(evaluation_point);
    return;
  }
};

// This specialization is needed to get a specialized recursion anchor for the
// case P = 0.
template <int P>
class HiddenBernsteinClass<-1, P> {
 public:
  static inline void EvalBasis(double *in,
                               const double evaluation_point) noexcept {
    BEMBEL_UNUSED_(in);
    BEMBEL_UNUSED_(evaluation_point);
    assert(
        false &&
        "Pos.C This should not happen. Something is wrong with the recursion");
  };
  static inline void EvalDerBasis(double *in,
                                  const double evaluation_point) noexcept {
    BEMBEL_UNUSED_(in);
    BEMBEL_UNUSED_(evaluation_point);
    assert(
        false &&
        "Pos.C This should not happen. Something is wrong with the recursion");
  };
};
////////////////////////////////////////////////////////////////////////////////
/// Evaluation Routines
////////////////////////////////////////////////////////////////////////////////

// think twice when replacing double* by something else. All other attempts
// were slower
template <int P>
void EvalBernsteinBasis(double *in, const double evaluation_point) noexcept {
  HiddenBernsteinClass<P, P>::EvalBasis(in, evaluation_point);
}

// think twice when replacing double* by something else. All other attempts
// were slower
template <int P>
void EvalBernsteinDerBasis(double *in, const double evaluation_point) noexcept {
  in[P] = P * Bernstein<P - 1, P - 1>(evaluation_point);
  HiddenBernsteinClass<P - 1, P - 1>::EvalDerBasis(in, evaluation_point);
  return;
}

}  // namespace Basis
}  // namespace Bembel
#endif  // BEMBEL_SRC_SPLINE_BERNSTEIN_HPP_
