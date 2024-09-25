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
 *  \brief computes a ball enclosing the union of \f$B_r1(mp1)\f$ and
 * \f$B_r2(mp2)\f$, i.e \f$B(mp,r)\supset B_r1(mp1) \cup B_r2(mp2)\f$.
 */
void computeEnclosingBall(Eigen::Vector3d *mp, double *r,
                          const Eigen::Vector3d &mp1, double r1,
                          const Eigen::Vector3d &mp2, double r2) {
  // compute distance vector of the two spheres
  auto z = mp1 - mp2;
  auto norm = (mp1 - mp2).norm();
  // B(d2,r2) subset B(d1,r1)
  if (norm + r2 <= r1) {
    *mp = mp1;
    *r = r1;
    // B(d1,r1) subset B(d2,r2)
  } else if (norm + r1 <= r2) {
    *mp = mp2;
    *r = r2;
    // the union is not a ball
  } else {
    *mp = 0.5 * (mp1 + mp2 + (r1 - r2) / norm * z);
    *r = 0.5 * (r1 + r2 + norm);
    *r = 0.5 * (r1 + r2 + norm);
  }
  return;
}
}  // namespace util
}  // namespace Bembel
#endif  // BEMBEL_SRC_UTIL_GEOMETRYHELPER_HPP_
