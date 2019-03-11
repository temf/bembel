// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef _shredder_
#define _shredder_

#include "patch.h"
namespace Spl {

/**
 *  @brief This is patch-shredder. It eats a single patch, chews it carefully,
 * and spits out a list of patches which are all in Bezier-Form.
 *
 */
inline std::vector<Patch> PatchShredder(const Patch &in) noexcept {
  // std::cout << "This is PatchShredder. I'm in." << std::endl;

  if (in.yuniqueknt.size() == 2 && in.xuniqueknt.size() == 2) {
    // std::cout << "Nothing to shred.\n";
    return {in};
  }

  const int xchips = in.xuniqueknt.size() - 1;
  const int ychips = in.yuniqueknt.size() - 1;

  const int xp = in.xp;
  const int yp = in.yp;
  // const bool turnNormal = in.turnNormal;
  // const bool flip = in.flip;
  const int numy = ychips * yp;

  std::vector<Patch> out(xchips * ychips);

  // std::cout << "What is going on...?" << std::endl;

  for (int ix = 0; ix < xchips; ix++) {
    for (int iy = 0; iy < ychips; iy++) {
      const int idx = ix * ychips + iy;

      // std::cout << "This is my " << idx << "st go at this" <<
      // std::endl;
      // Steal Settings of input
      out[idx].xuniqueknt = {0, 1};
      out[idx].yuniqueknt = {0, 1};
      out[idx].xp = xp;
      out[idx].yp = yp;
      // out[idx].turnNormal = turnNormal;
      // out[idx].flip = flip;
      out[idx].data.reserve(xp * yp * 4);
    }
  }

  // std::cout << "Initialized the stuff." << std::endl;

  for (int ix = 0; ix < xchips; ix++) {
    for (int iy = 0; iy < ychips; iy++) {
      const int idx = ix * ychips + iy;

      // This takes (x,y,z,w) * xp * yp control points for the
      // desired Bezier-Patch

      for (int jx = 0; jx < xp; jx++) {
        for (int jy = 0; jy < yp; jy++) {
          const int accs = 4 * (numy * (xp * ix + jx) + yp * iy + jy);
          for (int k = 0; k < 4; k++) {
            out[idx].data.push_back(in.data[accs + k]);
          }
        }
      }
    }
  }

  return out;
}

// Shredds a whole vector of Patches, oh my....
inline std::vector<Patch> PatchShredder(const std::vector<Patch> &ins) {
  std::vector<Patch> out;
  const int insize = ins.size();

  for (int i = 0; i < insize; i++) {
    std::vector<Patch> tmp = PatchShredder(ins[i]);

    const int tmpsize = tmp.size();

    for (int j = 0; j < tmpsize; j++) out.push_back(tmp[j]);
  }

  out.shrink_to_fit();

  return out;
}
}  // namespace Spl
#endif
