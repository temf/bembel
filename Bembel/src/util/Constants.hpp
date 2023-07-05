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
#ifndef BEMBEL_SRC_UTIL_CONSTANTS_HPP_
#define BEMBEL_SRC_UTIL_CONSTANTS_HPP_

namespace Bembel {
namespace Constants {

////////////////////////////////////////////////////////////////////////////////
/// variables
////////////////////////////////////////////////////////////////////////////////
constexpr double generic_tolerance = 1e-6;
// some not further specified constant
constexpr int MaxP = 20;

// some not further specified constant
constexpr int maximum_quadrature_degree = 50;
// constants for the mesh refinement
constexpr double corners[2][4] = {{0., 1., 1., 0}, {0., 0., 1., 1.}};
constexpr double llcs[2][4] = {{0., .5, .5, .0}, {0., 0., .5, .5}};
constexpr double edgemps[2][5] = {{.5, 1., .5, 0, .5}, {0., .5, 1., .5, .5}};
// tolerance for point comparison to determine patch topology
constexpr double pt_comp_tolerance = 1e-9;
// realloc size must be bigger than 4
constexpr size_t Bembel_alloc_size = 100;
// the interpolation problem solved during the assembly of the projector needs
// to filter some almost-zero coefficients that might be introduced during the
// solution of the linear system
constexpr double projector_tolerance = 1e-4;
////////////////////////////////////////////////////////////////////////////////
/// methods
////////////////////////////////////////////////////////////////////////////////

inline constexpr bool isAlmostZero(double in) {
  return in < generic_tolerance && in > -generic_tolerance;
}
}  // namespace Constants
}  // namespace Bembel

#endif  // BEMBEL_SRC_UTIL_CONSTANTS_HPP_
