// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2024 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_SRC_UTIL_GEOMETRYHELPER_HPP_
#define BEMBEL_SRC_UTIL_GEOMETRYHELPER_HPP_

namespace Bembel {
namespace util {

/**
 *  \brief computes a ball enclosing of two given enclosing balls.
 *
 * The computed enclosing ball contains for sure the two given balls with
 * mit_point1 and radius1 and 2, respectively.
 *
 * \param *midpoint: Pointer to Eigen::Vector3d which is the new mid point of
 * the enclosing ball
 * \param *radius: Pointer where the radius of new enclosing ball is saved.
 * \param midpoint1: Eigen::Vector3d midpoint of the first ball
 * \param radius1: double radius of the first ball
 * \param midpoint2: Eigen::Vector3d midpoint of the second ball
 * \param radius2: double radius of the second ball
 */
void computeEnclosingBall(Eigen::Vector3d *midpoint, double *radius,
                          const Eigen::Vector3d &midpoint1, double radius1,
                          const Eigen::Vector3d &midpoint2, double radius2) {
  // compute distance vector of the two spheres
  auto z = midpoint1 - midpoint2;
  auto norm = (midpoint1 - midpoint2).norm();
  // B(d2,radius2) subset B(d1,radius1)
  if (norm + radius2 <= radius1) {
    *midpoint = midpoint1;
    *radius = radius1;
    // B(d1,radius1) subset B(d2,radius2)
  } else if (norm + radius1 <= radius2) {
    *midpoint = midpoint2;
    *radius = radius2;
    // the union is not a ball
  } else {
    *midpoint = 0.5 * (midpoint1 + midpoint2 + (radius1 - radius2) / norm * z);
    *radius = 0.5 * (radius1 + radius2 + norm);
  }
  return;
}
}  // namespace util
}  // namespace Bembel
#endif  // BEMBEL_SRC_UTIL_GEOMETRYHELPER_HPP_
