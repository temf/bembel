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

#include <Bembel/Geometry>
#include <iostream>

int main() {
  Bembel::Geometry geometry("sphere.dat");
  std::cout << "The geometry has " << geometry.get_number_of_patches()
            << " patches." << std::endl;
  std::cout << std::string(60, '-') << std::endl;
  return 0;
}
