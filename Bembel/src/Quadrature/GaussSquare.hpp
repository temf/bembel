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

#ifndef BEMBEL_SRC_QUADRATURE_GAUSSSQUARE_HPP_
#define BEMBEL_SRC_QUADRATURE_GAUSSSQUARE_HPP_
namespace Bembel {

template <unsigned int Order>
using GaussSquare = TensorProductQuadratureVector<GaussLegendreRule, Order>;

}
#endif  // BEMBEL_SRC_QUADRATURE_GAUSSSQUARE_HPP_
