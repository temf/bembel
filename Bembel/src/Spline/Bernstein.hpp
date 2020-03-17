// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_SPLINE_BERNSTEIN_H_
#define BEMBEL_SPLINE_BERNSTEIN_H_

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
template <typename T, int N, int P>
class HiddenBernsteinClass {
 public:
  static inline T EvalCoefs(T *in, double evaluation_point) noexcept {
    return in[N] * Bernstein<N, P>(evaluation_point) +
           HiddenBernsteinClass<T, N - 1, P>::EvalCoefs(in, evaluation_point);
  }
  static inline T EvalDerCoefs(T *in, double evaluation_point) noexcept {
    return (
        (in[N + 1] - in[N]) * Bernstein<N, P>(evaluation_point) +
        HiddenBernsteinClass<T, N - 1, P>::EvalDerCoefs(in, evaluation_point));
  }
  static inline void EvalBasisPEQ(T *in, double evaluation_point) noexcept {
    in[N] += Bernstein<N, P>(evaluation_point);
    HiddenBernsteinClass<T, N - 1, P>::EvalBasisPEQ(in, evaluation_point);
    return;
  }
  static inline void EvalDerBasisPEQ(T *in, double evaluation_point) noexcept {
    in[N] += (P + 1) * (Bernstein<N - 1, P>(evaluation_point) -
                        Bernstein<N, P>(evaluation_point));
    HiddenBernsteinClass<T, N - 1, P>::EvalDerBasisPEQ(in, evaluation_point);
    return;
  }
  static inline void EvalBasis(T *in, double evaluation_point) noexcept {
    in[N] = Bernstein<N, P>(evaluation_point);
    HiddenBernsteinClass<T, N - 1, P>::EvalBasis(in, evaluation_point);
    return;
  }
  static inline void EvalDerBasis(T *in, double evaluation_point) noexcept {
    in[N] = (P + 1) * (Bernstein<N - 1, P>(evaluation_point) -
                       Bernstein<N, P>(evaluation_point));
    HiddenBernsteinClass<T, N - 1, P>::EvalDerBasis(in, evaluation_point);
    return;
  }
};

template <typename T, int P>
class HiddenBernsteinClass<T, 0, P> {
 public:
  static inline T EvalCoefs(T *in, double evaluation_point) noexcept {
    return in[0] * Bernstein<0, P>(evaluation_point);
  }
  static inline T EvalDerCoefs(T *in, double evaluation_point) noexcept {
    // P needs to be passed lower to avoid infinite recursion
    return (in[1] - in[0]) * Bernstein<0, P>(evaluation_point);
  }
  static inline void EvalBasisPEQ(T *in, double evaluation_point) noexcept {
    in[0] += Bernstein<0, P>(evaluation_point);
    return;
  }
  static inline void EvalDerBasisPEQ(T *in, double evaluation_point) noexcept {
    // P needs to be passed lower to avoid infinite recursion
    in[0] += (-P - 1) * Bernstein<0, P>(evaluation_point);
    return;
  }

  static inline void EvalBasis(T *in, double evaluation_point) noexcept {
    in[0] = Bernstein<0, P>(evaluation_point);
    return;
  }
  static inline void EvalDerBasis(T *in, double evaluation_point) noexcept {
    // P needs to be passed lower to avoid infinite recursion
    in[0] = (-P - 1) * Bernstein<0, P>(evaluation_point);
    return;
  }
};

// This specialization is needed to get a specialized recursion anchor for the
// case P = 0.
template <typename T, int P>
class HiddenBernsteinClass<T, -1, P> {
 public:
  static inline T EvalCoefs(T *in, double evaluation_point) noexcept {
    (void)in;
    (void)evaluation_point;
    assert(
        false &&
        "Pos.A This should not happen. Something is wrong with the recursion");
  };
  static inline T EvalDerCoefs(T *in, double evaluation_point) noexcept {
    // P needs to be passed lower to avoid infinite recursion
    (void)in;
    (void)evaluation_point;
    return 0;
  };
  static inline void EvalBasis(T *in, double evaluation_point) noexcept {
    (void)in;
    (void)evaluation_point;
    assert(
        false &&
        "Pos.C This should not happen. Something is wrong with the recursion");
  };
  static inline void EvalDerBasis(T *in, double evaluation_point) noexcept {
    (void)in;
    (void)evaluation_point;
    // P needs to be passed lower to avoid infinite recursion
    return;
  };
  static inline void EvalBasisPEQ(T *in, double evaluation_point) noexcept {
    (void)in;
    (void)evaluation_point;
    assert(
        false &&
        "Pos.C This should not happen. Something is wrong with the recursion");
  };
  static inline void EvalDerBasisPEQ(T *in, double evaluation_point) noexcept {
    (void)in;
    (void)evaluation_point;
    // P needs to be passed lower to avoid infinite recursion
    return;
  };
};
////////////////////////////////////////////////////////////////////////////////
/// Evaluation Routines
////////////////////////////////////////////////////////////////////////////////
template <typename T, int P>
T EvalBernstein(T *in, double evaluation_point) noexcept {
  return HiddenBernsteinClass<T, P, P>::EvalCoefs(in, evaluation_point);
}

template <typename T, int P>
void EvalBernstein(T *in, const std::vector<double> &evaluation_points,
                   T *out) noexcept {
  const int N = evaluation_points.size();
  for (int i = 0; i < N; i++)
    out[i] = HiddenBernsteinClass<T, P, P>::EvalCoefs(in, evaluation_points[i]);
  return;
}

template <typename T, int P>
std::vector<T> EvalBernstein(
    T *in, const std::vector<double> &evaluation_points) noexcept {
  const int N = evaluation_points.size();
  std::vector<double> out(N);
  for (int i = 0; i < N; i++)
    out[i] = HiddenBernsteinClass<T, P, P>::EvalCoefs(in, evaluation_points[i]);
  return out;
}

template <typename T, int P>
void EvalBernsteinBasisPEQ(T *in, double evaluation_point) noexcept {
  HiddenBernsteinClass<T, P, P>::EvalBasisPEQ(in, evaluation_point);
  return;
}

template <typename T, int P>
void EvalBernsteinBasis(T *in, double evaluation_point) noexcept {
  HiddenBernsteinClass<T, P, P>::EvalBasis(in, evaluation_point);
  return;
}
////////////////////////////////////////////////////////////////////////////////
/// Evaluation of the Derivatives
////////////////////////////////////////////////////////////////////////////////
template <typename T, int P>
T EvalBernsteinDer(T *in, double evaluation_point) noexcept {
  return P * HiddenBernsteinClass<T, P - 1, P - 1>::EvalDerCoefs(
                 in, evaluation_point);
}

template <typename T, int P>
void EvalBernsteinDer(T *in, const std::vector<double> &evaluation_points,
                      T *out) noexcept {
  const int N = evaluation_points.size();
  for (int i = 0; i < N; i++)
    out[i] = P * HiddenBernsteinClass<T, P - 1, P - 1>::EvalDerCoefs(
                     in, evaluation_points[i]);
  return;
}

template <typename T, int P>
std::vector<T> EvalBernsteinDer(
    T *in, const std::vector<double> &evaluation_points) noexcept {
  const int N = evaluation_points.size();
  std::vector<double> out(N);
  for (int i = 0; i < N; i++)
    out[i] = P * HiddenBernsteinClass<T, P - 1, P - 1>::EvalDerCoefs(
                     in, evaluation_points[i]);
  return out;
}

template <typename T, int P>
void EvalBernsteinDerBasisPEQ(T *in, double evaluation_point) noexcept {
  in[P] += P * Bernstein<P - 1, P - 1>(evaluation_point);
  HiddenBernsteinClass<T, P - 1, P - 1>::EvalDerBasisPEQ(in, evaluation_point);
  return;
}

template <typename T, int P>
void EvalBernsteinDerBasis(T *in, double evaluation_point) noexcept {
  in[P] = P * Bernstein<P - 1, P - 1>(evaluation_point);
  HiddenBernsteinClass<T, P - 1, P - 1>::EvalDerBasis(in, evaluation_point);
  return;
}

}  // namespace Basis
}  // namespace Bembel
#endif
