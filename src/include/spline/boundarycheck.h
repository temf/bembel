// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef _INCLUDE_BOUNDARY_CHECK_
#define _INCLUDE_BOUNDARY_CHECK_

#include <array>
#include "spline/boundarymatch.h"
#include "spline/patch.h"

namespace Spl {

/*
    matchpn =>
    0 -> x = 0 kante
    1 -> y = 1 kante
    2 -> x = 1 kante
    3 -> y = 0 kante
*/
/**
 *  @brief This checks weather the geometry can be used for computations in our
 * code. Edges with different directions of parametrization are allowed, but the
 * orientation of every surface patch needs to be conforming, i.e., all normals
 * need to point in the right direction.
 *
 * If you want to understand what this function does, take a piece of paper and
 * do the following: Draw a 4x4 grid of two squares next to each other on the
 * paper. These each to correspond to two patches with common edge. If the
 * normal should point in the same direction, there are 16 possible ways to
 * lable the directions of parametrisations, where the common edge corresponds
 * to a pair one of the above matchpn conditions, i.e., which edge is the shared
 * edge on each patch. If we find a combination, it is enough to simply check
 * weather the direction of parametrisation along the edge is the same, an
 * information provided by the glue struct.
 *
 */
inline bool check_geometry(const std::vector<Patch> &geoms) {
  const auto gls = boundary_match(geoms);
  bool globalindicator = true;
  for (auto g : gls) {
    bool indicator = false;
    switch (g.matchp1) {
      case (0): {
        switch (g.matchp2) {
          case (0):
            indicator = g.reverse ? true : false;
            break;
          case (1):
            indicator = g.reverse ? true : false;
            break;
          case (2):
            indicator = g.reverse ? false : true;
            break;
          case (3):
            indicator = g.reverse ? false : true;
            break;
          default:
            break;
        }
      }; break;
      case (1): {
        switch (g.matchp2) {
          case (0):
            indicator = g.reverse ? true : false;
            break;
          case (1):
            indicator = g.reverse ? true : false;
            break;
          case (2):
            indicator = g.reverse ? false : true;
            break;
          case (3):
            indicator = g.reverse ? false : true;
            break;
          default:
            break;
        }
      }; break;
      case (2): {
        switch (g.matchp2) {
          case (0):
            indicator = g.reverse ? false : true;
            break;
          case (1):
            indicator = g.reverse ? false : true;
            break;
          case (2):
            indicator = g.reverse ? true : false;
            break;
          case (3):
            indicator = g.reverse ? true : false;
            break;
          default:
            break;
        }
      }; break;
      case (3): {
        switch (g.matchp2) {
          case (0):
            indicator = g.reverse ? false : true;
            break;
          case (1):
            indicator = g.reverse ? false : true;
            break;
          case (2):
            indicator = g.reverse ? true : false;
            break;
          case (3):
            indicator = g.reverse ? true : false;
            break;
          default:
            break;
        }
      }; break;
      default:
        break;
    }
    if(not(indicator)){
    	globalindicator = false;
    }
  }
  return globalindicator;
}
}  // namespace Spl
#endif