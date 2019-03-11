// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_GLUE__
#define __BEMBEL_GLUE__

#include <array>
#include <vector>
#include "Eigen/Sparse"
#include "Spline.hpp"
#include "discretization.h"

namespace Bembel {

/**
 * @brief This class takes care of identifying DOFs on different edges, which
 *must be identified with one another.
 **/
class Gluebucket {
 public:
  int numofpatches;
  int simpledofs;
  int splinedofsnonglue;
  std::vector<int> gluefrom;
  std::vector<int> glueto;
  std::vector<int> directions;
  // int maxdof;
  // std::vector<bool> takencareof;

  Gluebucket(discretization *ds, const std::vector<Spl::glue> &gls, int p,
             int M, int kntrep);

  int size() {
    assert(gluefrom.size() == directions.size() and
           directions.size() == glueto.size());
    return gluefrom.size();
  }

  // bool gonnaglue(int p, int i) {
  //   const int size = gluefrom.size();
  //   for (int j = 0; j < size; j++) {
  //     if (gluefrom[j] == p * i) {
  //       return true;
  //     }
  //   }
  //   return false;
  // }

  // It assumes some stuff on the problem, see the asserts in
  // make_glue_matrix_maxwell

  std::array<int, 2> getPartnerId(int i) {
    const int sz = size();
    for (int k = 0; k < sz; k++) {
      if (gluefrom[k] == i) {
        // if (glueto[k] > i) {
        return {glueto[k], directions[k]};
        // } else {
        // std::cout << "Whoopsydaisy!\n";
        // }
      }
    }
    for (int k = 0; k < sz; k++) {
      if (glueto[k] == i) {
        // if (gluefrom[k] > i) {
        return {gluefrom[k], directions[k]};
        // } else {
        // std::cout << "Whoopsydaisy!\n";
        // }
      }
    }
    return {-1, 0};
  }

  void checkTheBucket() {
    int maxdof = 0;
    for (auto x : gluefrom) {
      maxdof = maxdof > x ? maxdof : x;
    }
    for (auto x : glueto) {
      maxdof = maxdof > x ? maxdof : x;
    }

    std::vector<int> counter(maxdof + 1);
    // std::cout << "The largest dof is " << maxdof << "\n";
    for (int k = 0; k < size(); k++) {
      counter[gluefrom[k]] += 1;
      counter[glueto[k]] += 1;
    }

    for (int k = 0; k < maxdof + 1; k++) {
      if (counter[k] < 3) {
        // std::cout << "The problem is dof number " << k << "\n";
        assert(false && "Too many dofs glued at once");
      }
    }
  }

  void printTheBucket() {
    std::cout << "Pairs are \n";
    for (int k = 0; k < size(); k++) {
      std::cout << gluefrom[k] << " , " << glueto[k] << "\n";
    }
    return;
  }
};

std::vector<int> getGlueDofIds_Maxwell(const int largecomp, const int smallcomp,
                                       const int cse);

Eigen::SparseMatrix<double> make_glue_matrix_laplace(discretization *ds, int p,
                                                     int kntrep, int M,
                                                     const int maxsizeofmatrix);
Eigen::SparseMatrix<double> make_glue_matrix_maxwell(discretization *ds, int p,
                                                     int kntrep, int M,
                                                     const int maxsize);

}  // namespace Bembel
#endif
