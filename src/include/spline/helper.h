// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef _HELPER_INCLUDED_
#define _HELPER_INCLUDED_

#include <Eigen/Dense>
#include <vector>
namespace Spl {

/**
 *  @brief Tiny helper functions.
 */

template <typename T>
inline double mypow(T x, int i) {
  T out = x;
  for (int k = 0; k < i; k++) out *= x;
  return x;
}

// Unrolls. Y-Dir first.
template <typename T>
inline Eigen::Matrix<T, -1, 1> unroll(
    const Eigen::Matrix<T, -1, -1> &bigone) noexcept {
  const int nx = bigone.cols();
  const int ny = bigone.rows();
  // std::cout << "unroll\n";
  Eigen::Matrix<T, -1, 1> out(nx * ny, 1);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) out(i * ny + j) = bigone(j, i);
  }
  return out;
}
}  // namespace Spl
#endif
