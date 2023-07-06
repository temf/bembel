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
#ifndef BEMBEL_SRC_UTIL_MACROS_HPP_
#define BEMBEL_SRC_UTIL_MACROS_HPP_

#ifndef M_PI
#define BEMBEL_PI 3.14159265358979323846264338327950288
#else
#define BEMBEL_PI M_PI
#endif

#define BEMBEL_UNUSED_(x) (void)(x)

#define BEMBEL_SIGNUM_(x) ((x >= 0) ? 1 : -1)

#define BEMBEL_SQUARED_(x) ((x) * (x))

#endif  // BEMBEL_SRC_UTIL_MACROS_HPP_
