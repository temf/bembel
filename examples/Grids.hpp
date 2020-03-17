// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_GRIDS_
#define __BEMBEL_GRIDS_
#include <Eigen/Dense>
#include <tuple>
#include <vector>

namespace Bembel {
namespace Util {

inline Eigen::Matrix<double, Eigen::Dynamic, 3> makeTensorProductGrid(
    Eigen::VectorXd X, Eigen::VectorXd Y, Eigen::VectorXd Z) {
  const int maxX = std::max(X.rows(), X.cols());
  const int maxY = std::max(Y.rows(), Y.cols());
  const int maxZ = std::max(Z.rows(), Z.cols());
  Eigen::Matrix<double, Eigen::Dynamic, 3> out(maxX * maxY * maxZ, 3);

  for (int iz = 0; iz < maxZ; iz++)
    for (int iy = 0; iy < maxY; iy++)
      for (int ix = 0; ix < maxX; ix++) {
        out.row(ix + iy * maxX + iz * maxY * maxX) =
            Eigen::Vector3d(X(ix), Y(iy), Z(iz));
      }
  return out;
}

inline Eigen::Matrix<double, Eigen::Dynamic, 3> makeSphereGrid(
    const double r, const int n,
    const Eigen::Vector3d center = Eigen::Vector3d(0, 0, 0)) {
  Eigen::Matrix<double, Eigen::Dynamic, 3> out(n * n, 3);
  const double h = 1. / n;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      out.row(j + i * n) =
          (Eigen::Vector3d(
               r * cos(3.141592653 * h * i) * sin(3.141592653 * h * (j + 0.5)),
               r * sin(3.141592653 * h * i) * sin(3.141592653 * h * (j + 0.5)),
               r * cos(3.141592653 * h * j)) +
           center);
    }
  }

  return out;
}
}  // namespace Util
}  // namespace Bembel
#endif
