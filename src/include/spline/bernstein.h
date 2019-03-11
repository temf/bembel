// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef _bernstein_included_
#define _bernstein_included_

#include <assert.h>
#include <array>
#include <cstdlib>
#include <vector>
#include "spline/pascal.h"
namespace Spl {
/**
 *  @brief Template recursion to produce Bernstein polynomials. This is only
 * limited by the binomial coefficient, see pascal.h
 *
 */

template <int N>
inline constexpr double _brnstnX(double in) {
#ifdef _spline_debug_flag_
  assert(in > -.0000001 and in < 1.0000001 and
         "Function only valid for 0 <= x <= 1!");
#endif
  return in * _brnstnX<N - 1>(in);
}

template <>
inline constexpr double _brnstnX<1>(double in) {
  return in;
}
template <>
inline constexpr double _brnstnX<0>(double in) {
  return 1.;
}
template <>
inline constexpr double _brnstnX<-1>(double in) {
  return 0.;
}

template <int N, int P>
inline constexpr double brnstn(double in) {
  return binomial<N, P>::value * _brnstnX<N>(in) * _brnstnX<P - N>(1. - in);
}

//// Hidden brnstn Classes

template <typename T, int N, int P>
class _HiddenbrnstnCls {
 public:
  static inline T evalCoefs(T *in, double x) {
    return in[N] * brnstn<N, P>(x) +
           _HiddenbrnstnCls<T, N - 1, P>::evalCoefs(in, x);
  };
  static inline T evalDerCoefs(T *in, double x) {
    return ((in[N + 1] - in[N]) * brnstn<N, P>(x) +
            _HiddenbrnstnCls<T, N - 1, P>::evalDerCoefs(in, x));
  };
  static inline void evalBasisPEQ(T *in, double x) {
    in[N] += brnstn<N, P>(x);
    _HiddenbrnstnCls<T, N - 1, P>::evalBasisPEQ(in, x);
    return;
  };
  static inline void evalDerBasisPEQ(T *in, double x) {
    in[N] += (P + 1) * (brnstn<N - 1, P>(x) - brnstn<N, P>(x));
    _HiddenbrnstnCls<T, N - 1, P>::evalDerBasisPEQ(in, x);
    return;
  };
  static inline void evalBasis(T *in, double x) {
    in[N] = brnstn<N, P>(x);
    _HiddenbrnstnCls<T, N - 1, P>::evalBasis(in, x);
    return;
  };
  static inline void evalDerBasis(T *in, double x) {
    in[N] = (P + 1) * (brnstn<N - 1, P>(x) - brnstn<N, P>(x));
    _HiddenbrnstnCls<T, N - 1, P>::evalDerBasis(in, x);
    return;
  };
};

template <typename T, int P>
class _HiddenbrnstnCls<T, 0, P> {
 public:
  static inline T evalCoefs(T *in, double x) {
    return in[0] * brnstn<0, P>(x);
  };
  static inline T evalDerCoefs(T *in, double x) {
    // P needs to be passed lower to avoid infinite recursion
    return (in[1] - in[0]) * brnstn<0, P>(x);
  };
  static inline void evalBasisPEQ(T *in, double x) {
    in[0] += brnstn<0, P>(x);
    return;
  };
  static inline void evalDerBasisPEQ(T *in, double x) {
    // P needs to be passed lower to avoid infinite recursion
    in[0] += (-P - 1) * brnstn<0, P>(x);
    return;
  };

  static inline void evalBasis(T *in, double x) {
    in[0] = brnstn<0, P>(x);
    return;
  };
  static inline void evalDerBasis(T *in, double x) {
    // P needs to be passed lower to avoid infinite recursion
    in[0] = (-P - 1) * brnstn<0, P>(x);
    return;
  };
};

// This specialization is needed to get a specialized recursion anchor for the
// case P = 0.
template <typename T, int P>
class _HiddenbrnstnCls<T, -1, P> {
 public:
  static inline T evalCoefs(T *in, double x) {
    (void)in;
    (void)x;
    assert(
        false &&
        "Pos.A This should not happen. Something is wrong with the recursion");
  };
  static inline T evalDerCoefs(T *in, double x) {
    // P needs to be passed lower to avoid infinite recursion
    (void)in;
    (void)x;
    return 0;
  };
  static inline void evalBasis(T *in, double x) {
    (void)in;
    (void)x;
    assert(
        false &&
        "Pos.C This should not happen. Something is wrong with the recursion");
  };
  static inline void evalDerBasis(T *in, double x) {
    (void)in;
    (void)x;
    // P needs to be passed lower to avoid infinite recursion
    return;
  };
  static inline void evalBasisPEQ(T *in, double x) {
    (void)in;
    (void)x;
    assert(
        false &&
        "Pos.C This should not happen. Something is wrong with the recursion");
  };
  static inline void evalDerBasisPEQ(T *in, double x) {
    (void)in;
    (void)x;
    // P needs to be passed lower to avoid infinite recursion
    return;
  };
};

////// Evaluation Routines

template <typename T, int P>
T evalBrnstn(T *in, double x) {
  return _HiddenbrnstnCls<T, P, P>::evalCoefs(in, x);
}

template <typename T, int P>
void evalBrnstn(T *in, const std::vector<double> &x, T *out) {
  const int N = x.size();
  for (int i = 0; i < N; i++)
    out[i] = _HiddenbrnstnCls<T, P, P>::evalCoefs(in, x[i]);
  return;
}

template <typename T, int P>
std::vector<T> evalBrnstn(T *in, const std::vector<double> &x) {
  const int N = x.size();
  std::vector<double> out(N);
  for (int i = 0; i < N; i++)
    out[i] = _HiddenbrnstnCls<T, P, P>::evalCoefs(in, x[i]);
  return out;
}

template <typename T, int P>
void evalBrnstnBasisPEQ(T *in, double x) {
  _HiddenbrnstnCls<T, P, P>::evalBasisPEQ(in, x);
  return;
}

template <typename T, int P>
void evalBrnstnBasis(T *in, double x) {
  _HiddenbrnstnCls<T, P, P>::evalBasis(in, x);
  return;
}

///////// Evaluation of the Derivatives

template <typename T, int P>
T evalBrnstnDer(T *in, double x) {
  return P * _HiddenbrnstnCls<T, P - 1, P - 1>::evalDerCoefs(in, x);
}

template <typename T, int P>
void evalBrnstnDer(T *in, const std::vector<double> &x, T *out) {
  const int N = x.size();
  for (int i = 0; i < N; i++)
    out[i] = P * _HiddenbrnstnCls<T, P - 1, P - 1>::evalDerCoefs(in, x[i]);
  return;
}

template <typename T, int P>
std::vector<T> evalBrnstnDer(T *in, const std::vector<double> &x) {
  const int N = x.size();
  std::vector<double> out(N);
  for (int i = 0; i < N; i++)
    out[i] = P * _HiddenbrnstnCls<T, P - 1, P - 1>::evalDerCoefs(in, x[i]);
  return out;
}

template <typename T, int P>
void evalBrnstnDerBasisPEQ(T *in, double x) {
  in[P] += P * brnstn<P - 1, P - 1>(x);
  _HiddenbrnstnCls<T, P - 1, P - 1>::evalDerBasisPEQ(in, x);
  return;
}

template <typename T, int P>
void evalBrnstnDerBasis(T *in, double x) {
  in[P] = P * brnstn<P - 1, P - 1>(x);
  _HiddenbrnstnCls<T, P - 1, P - 1>::evalDerBasis(in, x);
  return;
}

}  // namespace Spl
#endif