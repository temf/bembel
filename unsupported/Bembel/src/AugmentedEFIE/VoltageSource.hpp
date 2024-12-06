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
#ifndef BEMBEL_SRC_AUGMENTEDEFIE_VOLTAGESOURCE_HPP_
#define BEMBEL_SRC_AUGMENTEDEFIE_VOLTAGESOURCE_HPP_

namespace Bembel {
/**
 *  \brief Struct to handle the geometry information of the excitation by a
 *  voltage source.
 *
 *  The vector contains all patches of the positive and negative pole
 *  respectively. The first entry of the array is the patch number whereas the
 *  second entry corresponds to the edge case which is connected to the voltage
 *  source.
 */
struct VoltageSource {
  std::vector<std::array<int, 2>> positive_;
  std::vector<std::array<int, 2>> negative_;
};
}  // namespace Bembel
#endif  // BEMBEL_SRC_AUGMENTEDEFIE_VOLTAGESOURCE_HPP_
