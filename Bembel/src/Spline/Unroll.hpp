// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_SPLINE_HELPER_H_
#define BEMBEL_SPLINE_HELPER_H_

namespace Bembel {
namespace Spl {

/**
 *  \brief Tiny helper functions.
 */

// Unrolls a matrix into a vector. Y-Dir first.
template <typename T>
inline Eigen::Matrix<T, -1, 1> Unroll(
    const Eigen::Matrix<T, -1, -1> &input_matrix) noexcept {
  const int nx = input_matrix.cols();
  const int ny = input_matrix.rows();
  // std::cout << "unroll\n";
  Eigen::Matrix<T, -1, 1> out(nx * ny, 1);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) out(i * ny + j) = input_matrix(j, i);
  }
  return out;
}
}  // namespace Spl
}  // namespace Bembel
#endif  // BEMBEL_SPLINE_HELPER_H_
