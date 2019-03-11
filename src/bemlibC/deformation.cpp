// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "Geometry.hpp"
#include "Spline.hpp"

namespace Bembel {

// Functions to deform the geometry.
void Geometry::rotate(double x, double y, double z) {
  Eigen::Matrix<double, 3, 3> trafo = Eigen::Matrix<double, 3, 3>::Zero(3, 3);
  trafo(0, 0) = 1;
  trafo(1, 1) = cos(x);
  trafo(2, 2) = cos(x);
  trafo(1, 2) = -sin(x);
  trafo(2, 1) = sin(x);

  // y-axis rotation
  Eigen::Matrix<double, 3, 3> rotmat = Eigen::Matrix<double, 3, 3>::Zero(3, 3);
  rotmat(0, 0) = cos(y);
  rotmat(1, 1) = 1;
  rotmat(2, 2) = cos(y);
  rotmat(0, 2) = -sin(y);
  rotmat(2, 0) = sin(y);
  trafo = (rotmat * trafo).eval();

  // z-axis rotation
  rotmat = Eigen::Matrix<double, 3, 3>::Zero(3, 3);
  rotmat(0, 0) = cos(z);
  rotmat(1, 1) = cos(z);
  rotmat(0, 1) = -sin(z);
  rotmat(1, 0) = sin(z);
  rotmat(2, 2) = 1;
  trafo = (rotmat * trafo).eval();

  deform(trafo);

  return;
}

void Geometry::stretch(double x, double y, double z) {
  Eigen::Matrix<double, 3, 3> stretchmat =
      Eigen::Matrix<double, 3, 3>::Zero(3, 3);

  stretchmat(0, 0) = x;
  stretchmat(1, 1) = y;
  stretchmat(2, 2) = z;

  deform(stretchmat);

  return;
}

void Geometry::deform(Eigen::Matrix<double, 3, 3> trafo) {
  const int num_of_patches = _geom.size();
  for (int i = 0; i < num_of_patches; i++) {
    const int size = _geom[i].data.size();
    const int pts = size / 4;
    for (int j = 0; j < pts; j++) {
      // Note here that k<3 instead of k<4, since the
      // case k==3 is covered because the weight does not change
      Eigen::Vector3d vec;
      for (int k = 0; k < 3; k++) {
        vec(k) = _geom[i].data[j * 4 + k];
      }

      vec = (trafo * vec).eval();

      for (int k = 0; k < 3; k++) {
        _geom[i].data[j * 4 + k] = vec(k);
      }
    }
  }
  return;
}
}  // namespace Bembel
