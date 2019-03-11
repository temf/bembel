// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "glue.h"

namespace Bembel {
Eigen::SparseMatrix<double> make_glue_matrix_laplace(
    discretization *ds, int p, int kntrep, int M, const int maxsizeofmatrix) {
  auto gls = boundary_match(*ds->mesh->geom);

  const int a_o = ds->a_o;
  const int dofsX = (a_o) + ((1 << M) - 1) * std::min(a_o, kntrep);
  const int dofsY = (a_o) + ((1 << M) - 1) * std::min(a_o, kntrep);
  const int dofs = dofsX * dofsY;

  std::vector<int> gluefrom;
  std::vector<int> glueto;
  std::vector<int> from;
  std::vector<int> to;
  from.reserve(dofsX);
  to.reserve(dofsX);

  auto getDofId = [&](int j, int matchp) {
    // Now, we look at the dof numbers of the individual patch.
    // This switch does some arithmetic: The cases are :
    // 0 : x == 0 Kante
    // 1 : Y == 1 Kante
    // 2 : X == 1 Kante
    // 3 : Y == 0 Kante
    // See boundarymatch.h for more info.

    switch (matchp) {
      case (0):
        return j * dofsX;
      case (1):
        return dofs - dofsX + j;
      case (2):
        return (j + 1) * dofsX - 1;
      default:
        return j;
    }
  };

  for (auto g : gls) {
    for (int j = 0; j < dofsX; j++) {
      from.push_back(getDofId(j, g.matchp1) + g.patch1 * dofs);
      to.push_back(getDofId(j, g.matchp2) + g.patch2 * dofs);
    }

    for (int i = 0; i < dofsX; i++) {
      gluefrom.push_back(from[i]);
      glueto.push_back(to[(g.reverse ? (dofsX - i - 1) : i)]);
    }
  }

  assert(glueto.size() == gluefrom.size());

  auto getPartnerId = [&](int i) {
    const int sz = gluefrom.size();
    for (int k = 0; k < sz; k++) {
      if (gluefrom[k] == i) {
        if (glueto[k] > i) {
          return glueto[k];
        } else {
        }
      }
    }
    for (int k = 0; k < sz; k++) {
      if (glueto[k] == i) {
        if (gluefrom[k] > i) {
          return gluefrom[k];
        } else {
        }
      }
    }
    return -1;
  };

  std::vector<int> takencareof;
  std::vector<Eigen::Triplet<int>> trips;

  auto isTakenCareOf = [&takencareof](int in) {
    for (auto x : takencareof) {
      if (x == in) {
        return true;
      }
    }
    return false;
  };

  int count = 0;
  for (int i = 0; i < maxsizeofmatrix; i++) {
    if (isTakenCareOf(i)) {
      // std::cout << "I skip  i = " << i << "\n";
      count = count + 1;
    } else {
      // std::cout << "I identify (" << i - count << "," << i << "," << 1 <<
      // ")\n";
      trips.push_back(Eigen::Triplet<int>(i - count, i, 1));
      // This is negative if there is no partner!
      const int pid = getPartnerId(i);
      if (pid >= 0) {
        takencareof.push_back(pid);
        trips.push_back(Eigen::Triplet<int>(i - count, pid, 1));
        // std::cout << "I glue     (" << i - count << "," << pid[0] << ","
        // << pid[1] << ")\n";
      }
    }
  }

  Eigen::SparseMatrix<double> gluemat(maxsizeofmatrix - count, maxsizeofmatrix);
  gluemat.setFromTriplets(trips.begin(), trips.end());

  return gluemat;
}

}  // namespace Bembel