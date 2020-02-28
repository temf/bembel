// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_SPLINE_PASCAL_H_
#define BEMBEL_SPLINE_PASCAL_H_

namespace Bembel {

/**
 *  \ingroup Spline
 *  \brief Here we hardcode the binomial coefficients through template
 *         recursion and perform bounds checking, just in case...
 *
 */
template <bool condition>
struct static_assertKlessseqN {};
template <>
struct static_assertKlessseqN<true> {
  enum { YOU_CALLED_BINOMIAL_WITH_KgtrN_OR_NEGATIVE_VALUES = 1 };
};

template <int K, int N>
struct Binomial {
  enum {
    value = static_assertKlessseqN<(K >= 0 && N >= 0 && K <= N)>::
                    YOU_CALLED_BINOMIAL_WITH_KgtrN_OR_NEGATIVE_VALUES
                ? Binomial<K - 1, N - 1>::value + Binomial<K, N - 1>::value
                : 0
  };
};
template <>
struct Binomial<0, 0> {
  enum { value = 1 };
};
template <int N>
struct Binomial<N, N> {
  enum { value = 1 };
};
template <int N>
struct Binomial<0, N> {
  enum { value = 1 };
};

}  // namespace Bembel
#endif
