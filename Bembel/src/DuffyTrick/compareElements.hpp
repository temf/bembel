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

#ifndef BEMBEL_SRC_DUFFYTRICK_COMPAREELEMENTS_HPP_
#define BEMBEL_SRC_DUFFYTRICK_COMPAREELEMENTS_HPP_

namespace Bembel {
namespace DuffyTrick {
/**
 *  \ingroup DuffyTrick
 *  \brief compares two elements for similarities and determines, how the
 *         elements have to be rotated to move the similarity to the first
 *         vertices_ or edge
 **/
Eigen::Vector3i compareElements(const ElementTreeNode &e1,
                                const ElementTreeNode &e2, double *dist) {
  Eigen::Vector3i retval;
  retval.setZero();
  // check if the two elements are identical and directly return;
  if (std::addressof(e1) == std::addressof(e2)) {
    retval << 4, 4, 2;
    *dist = 0;
    return retval;
  } else {
    // if they are not identical, check if they have a positive distance
    // check for common vertices
    *dist = (e1.midpoint_ - e2.midpoint_).norm() - e1.radius_ - e2.radius_;
    *dist = *dist >= 0 ? *dist : 0;
    // check if elements are distinct and return now
    if (*dist > .5 / (1 << e1.level_)) {
      retval << 4, 4, 0;
      return retval;
      // otherwise check for common edge/vertex. Note that there is a
      // short circuit: either two elements share a single edge or
      // single point. everything else will break the code
    } else {
      for (auto rot1 = 0; rot1 < 4; ++rot1)
        for (auto rot2 = 0; rot2 < 4; ++rot2)
          // check for common vertices_
          if (e1.vertices_[rot1] == e2.vertices_[rot2]) {
            // if there is a common vertices_, check for common edge
            if (e1.vertices_[3] == e2.vertices_[(rot2 + 1) % 4]) {
              retval << 3, rot2, 3;
              return retval;
            } else if (e1.vertices_[(rot1 + 1) % 4] ==
                       e2.vertices_[(rot2 + 3) % 4]) {
              retval << rot1, (rot2 + 3) % 4, 3;
              return retval;
            } else {
              retval << rot1, rot2, 4;
              return retval;
            }
          }
      retval << 4, 4, 0;
      return retval;
    }
  }
  // if you ended up here, something went terribly wrong
  retval << 4, 4, -1;
  return retval;
}
}  // namespace DuffyTrick
}  // namespace Bembel

#endif  // BEMBEL_SRC_DUFFYTRICK_COMPAREELEMENTS_HPP_
