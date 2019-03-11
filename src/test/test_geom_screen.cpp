// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "bemtest.h"

namespace Bembel {
namespace Test {

std::vector<Spl::Patch> mkScreen() {
  Spl::Patch bot;

  std::vector<double> tempknt = {0, 0, 1, 1};

  {
    Eigen::Matrix<double, 2, 2> X, Y, Z, ws;
    ws << 1, 1, 1, 1;
    X << 0, 1, 0, 1;
    Y << 0, 0, 0, 0;
    Z << 1, 1, 0, 0;

    X = X.rowwise().reverse().eval();
    Y = Y.rowwise().reverse().eval();
    Z = Z.rowwise().reverse().eval();
    ws = ws.rowwise().reverse().eval();
    std::vector<Eigen::Matrix<double, -1, -1>> tempctrl = {X, Y, Z, ws};
    bot.init(tempctrl, tempknt, tempknt);
  }

  std::vector<Spl::Patch> out = {bot};

  return out;
}

}  // namespace Test
}  // namespace Bembel