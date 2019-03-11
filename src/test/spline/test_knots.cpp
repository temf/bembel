// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "spltest.h"

int test_knots() {
  using namespace Spl;
  std::default_random_engine generator(time(0));
  std::uniform_int_distribution<int> distribution(1, 100);

  // Testing routines for knot-vec-generations and information extraction

  for (int i = 0; i < 200; i++) {
    int rand = distribution(generator);
    if (not(get_deg_of_knts(
                make_unif_knots(rand, rand + distribution(generator))) == rand))
      return 1;
  }

  for (int i = 1; i < 200; i++) {
    int rand = distribution(generator);
    if (not(extract_unique_knots(make_bezier_knots(distribution(generator)))
                .size() == 2))
      return 1;
    if (not(extract_unique_knots(make_unif_knots(distribution(generator), rand))
                .size() == 2 + rand))
      return 1;
  }

  return 0;
}