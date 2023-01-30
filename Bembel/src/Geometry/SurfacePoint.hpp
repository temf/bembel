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
//
#ifndef BEMBEL_SRC_GEOMETRY_SURFACEPOINT_HPP_
#define BEMBEL_SRC_GEOMETRY_SURFACEPOINT_HPP_
/**
 * \ingroup Geometry
 * \brief typedef of SurfacePoint
 *
 * This typedef is essential for any evaluation of a bilinear form. It provides
 * all required geometry information, stored as follows:
 * (0) x-coordinate of the evaluation point in the parameter domain [0,1]^2
 *     of the current element, i.e we map [0,1]^2->element->surface
 * (1) y-coordinate of the evaluation point in the parameter domain [0,1]^2
 *     of the current element, i.e we map [0,1]^2->element->surface
 * (2) a quadrature weight. Can be left empty if not used as part of a
 *     quadrature.
 * (3) x-coordinate of patch eval in space
 * (4) y-coordinate of patch eval in space
 * (5) z-coordinate of patch eval in space
 * (6) x-component of derivative in x-dir
 * (7) y-component of derivative in x-dir
 * (8) z-component of derivative in x-dir
 * (9) x-component of derivative in y-dir
 * (10) y-component of derivative in y-dir
 * (11) z-component of derivative in y-dir
 * For application of the pull-back to the reference domain, one requires the
 * jacobian of any point on the surface. Calling eval and evalJacobian of the
 * Patch class introduces work that needs to be done twice. The
 * updateSurdacePoint method is specialized and should be used, since it avoids
 * redundant work.
 **/
typedef Eigen::Matrix<double, 12, 1> SurfacePoint;

#endif  // BEMBEL_SRC_GEOMETRY_SURFACEPOINT_HPP_
