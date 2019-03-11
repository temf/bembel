// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef DSD_boundarymatch
#define DSD_boundarymatch

#include <array>
#include <vector>
#include "Eigen/Sparse"
#include "spline/patch.h"
/**
 *  @brief In this header file we have routines for the identifcation of common
 * edges and to give information about the topology to the gluebucket.
 *
 */
namespace Spl {
// This littel glue-struct stores all informations that is required about
// matching patches.

struct glue {
  int patch1;    //
  int patch2;    // Which two patches are neighbours?
  int gluecond;  // How is  the dir of param. wrt. each other?
  int matchp1;   //
  int matchp2;   // Which edges are the matching ones?
  bool reverse;
};

inline bool gluecheck(const glue &g) { return (g.patch1 < g.patch2); }

/*
    matchpn =>
    0 -> x = 0 kante
    1 -> y = 1 kante
    2 -> x = 1 kante
    3 -> y = 0 kante
*/

// The matching condition allows for nonmatching interfaces. This checks exactly
// that
inline bool need_reverse(glue g, const std::vector<Patch> &ps) {
  // int d[4] = {1, 0, 3, 2};
  // if ((g.matchp1 == g.matchp2) or (d[g.matchp1] == g.matchp2)) {
  //   // std::cout << "It happens! Patch " << g.patch1 << " and "<<g.patch2 <<
  //   // "\n"; assert(false);
  //   return true;
  // } else {
  //   return false;
  // }
  std::vector<double> p1p, p2p;
  switch (g.matchp1) {
    case (0):
      p1p = {0, 0, 0, 1};
      break;
    case (1):
      p1p = {0, 1, 1, 1};
      break;
    case (2):
      p1p = {1, 1, 0, 1};
      break;
    default:
      assert(g.matchp1 == 3);
      p1p = {0, 1, 0, 0};
      break;
  }

  switch (g.matchp2) {
    case (0):
      p2p = {0, 0, 0, 1};
      break;
    case (1):
      p2p = {0, 1, 1, 1};
      break;
    case (2):
      p2p = {1, 1, 0, 1};
      break;
    default:
      assert(g.matchp2 == 3);
      p2p = {0, 1, 0, 0};
      break;
  }

  Eigen::Vector3d x1 = ps[g.patch1].eval(p1p[0], p1p[2]);
  Eigen::Vector3d x2 = ps[g.patch1].eval(p1p[1], p1p[3]);

  Eigen::Vector3d y1 = ps[g.patch2].eval(p2p[0], p2p[2]);
  Eigen::Vector3d y2 = ps[g.patch2].eval(p2p[1], p2p[3]);

  constexpr double tol = 10e-9;

  if ((x1 - y1).norm() < tol and (x2 - y2).norm() < tol) {
    return false;
  } else {
    assert((x1 - y2).norm() < tol and (x2 - y1).norm() < tol);
    return true;
  }
}

inline std::vector<glue> boundary_match(const std::vector<Patch> &ps) {
  constexpr double tol = 10e-9;

  const int no_of_patches = ps.size();
  std::vector<glue> out;
  out.reserve(no_of_patches * 2);

  std::vector<std::array<Eigen::Matrix<double, 3, 1>, 4>> cornerinfo;
  cornerinfo.reserve(no_of_patches);

  // I precompute the corners of the ptches in accordance to matchpn, see
  // above.
  for (int i = 0; i < no_of_patches; i++) {
    std::array<Eigen::Matrix<double, 3, 1>, 4> pi;
    // pi[0] = (ps[i].eval(1.0, 0.0)); //
    // pi[1] = (ps[i].eval(0.0, 0.0)); //
    // pi[2] = (ps[i].eval(0.0, 1.0)); //
    // pi[3] = (ps[i].eval(1.0, 1.0)); //
    // pi[4] = (ps[i].eval(1.0, 0.0)); //
    // pi[5] = (ps[i].eval(0.0, 0.0)); //
    pi[0] = (ps[i].eval(0.0, 0.5));  // + ps[i].eval(0.0, 1.0));
    pi[1] = (ps[i].eval(0.5, 1.0));  // + ps[i].eval(1.0, 1.0));
    pi[2] = (ps[i].eval(1.0, 0.5));  // + ps[i].eval(1.0, 1.0));
    pi[3] = (ps[i].eval(0.5, 0.0));  // + ps[i].eval(1.0, 0.0));
    cornerinfo.push_back(pi);
  }

  // Debugging only
  // for (int i = 0; i < no_of_patches; i++) {
  // 	for (int j = 0; j < 5; j++)
  // 		std::cout << cornerinfo[i][j].transpose() << "\n";
  // 	std::cout << "\n";
  // }
  // end dbg

  for (int pi = 0; pi < no_of_patches; pi++) {
    // std::cout << "New pi \n\n";

    for (int pj = pi; pj < no_of_patches; pj++) {
      // std::cout << "   Loop " << pj << "\n";

      // 	break;

      auto compareedge = [&cornerinfo](int ppi, int ppj, int ii, int jj) {
        return ((cornerinfo[ppi][ii] - cornerinfo[ppj][jj]).norm() < tol);
      };

      auto setGlue = [&pi, &pj, &ps](int i, int j, int k) {
        glue tmp;
        tmp.gluecond = k;
        tmp.patch1 = pi;
        tmp.patch2 = pj;
        tmp.matchp1 = i;
        tmp.matchp2 = j;
        // This requires a geometry satisfying the matching condition to work
        tmp.reverse = need_reverse(tmp, ps);
        return tmp;
      };

      // In case one patch has an edge with itself, this loop is called
      if (pi == pj) {
#if 0
        bool searching = true;
        if (searching and
            ((cornerinfo[pi][0] - cornerinfo[pj][1]).norm() < tol) and
            ((cornerinfo[pi][2] - cornerinfo[pj][3]).norm() < tol)) {
          searching = false;
        }

        if (searching and
            ((cornerinfo[pi][1] - cornerinfo[pj][2]).norm() < tol) and
            ((cornerinfo[pi][3] - cornerinfo[pj][4]).norm() < tol)) {
          searching = false;
        }

        if (searching and
            ((cornerinfo[pi][2] - cornerinfo[pj][3]).norm() < tol) and
            ((cornerinfo[pi][0] - cornerinfo[pj][1]).norm() < tol)) {
          searching = false;
        }

        if (searching and
            ((cornerinfo[pi][3] - cornerinfo[pj][4]).norm() < tol) and
            ((cornerinfo[pi][1] - cornerinfo[pj][2]).norm() < tol)) {
          searching = false;
        }

        assert(searching == true &&
               "Attention: * * * * * * * * * * * * * * * *\n     A "
               "patch needs to be glued "
               "with itself.\n     boundarymatch.h is still "
               "missing glue\n     conditions for such a case!\n * "
               "* * * * * * * * * * * * * * * * * * * *\n");

#endif
      }
      // for other cases, we check here.
      else {
        bool searching = true;

        if (searching and compareedge(pi, pj, 0, 0)) {
          // std::cout << "   In some case0\n";
          out.push_back(setGlue(0, 0, -1));
          searching = false;
        }
        if (searching and compareedge(pi, pj, 0, 1)) {
          // std::cout << "   In some case0\n";
          out.push_back(setGlue(0, 1, 1));
          searching = false;
        }
        if (searching and compareedge(pi, pj, 0, 2)) {
          // std::cout << "   In some case0\n";
          out.push_back(setGlue(0, 2, 1));
          searching = false;
        }
        if (searching and compareedge(pi, pj, 0, 3)) {
          // std::cout << "   In some case0\n";
          out.push_back(setGlue(0, 3, -1));
          searching = false;
        }

        if (searching and compareedge(pi, pj, 1, 0)) {
          // std::cout << "   In some case\n";
          out.push_back(setGlue(1, 0, 1));
          searching = false;
        }
        if (searching and compareedge(pi, pj, 1, 1)) {
          // std::cout << "   In some case\n";
          out.push_back(setGlue(1, 1, -1));
          searching = false;
        }
        if (searching and compareedge(pi, pj, 1, 2)) {
          // std::cout << "   In some case\n";
          out.push_back(setGlue(1, 2, -1));
          searching = false;
        }
        if (searching and compareedge(pi, pj, 1, 3)) {
          // std::cout << "   In some case\n";
          out.push_back(setGlue(1, 3, 1));
          searching = false;
        }

        if (searching and compareedge(pi, pj, 2, 0)) {
          // std::cout << "   In some case\n";
          out.push_back(setGlue(2, 0, 1));
          searching = false;
        }
        if (searching and compareedge(pi, pj, 2, 1)) {
          // std::cout << "   In some case\n";
          out.push_back(setGlue(2, 1, -1));
          searching = false;
        }
        if (searching and compareedge(pi, pj, 2, 2)) {
          // std::cout << "   In some case\n";
          out.push_back(setGlue(2, 2, -1));
          searching = false;
        }
        if (searching and compareedge(pi, pj, 2, 3)) {
          // std::cout << "   In some case\n";
          out.push_back(setGlue(2, 3, 1));
          searching = false;
        }

        if (searching and compareedge(pi, pj, 3, 0)) {
          // std::cout << "   In some case3\n";
          out.push_back(setGlue(3, 0, -1));
          searching = false;
        }
        if (searching and compareedge(pi, pj, 3, 1)) {
          // std::cout << "   In some case3\n";
          out.push_back(setGlue(3, 1, 1));
          searching = false;
        }
        if (searching and compareedge(pi, pj, 3, 2)) {
          // std::cout << "   In some case3\n";
          out.push_back(setGlue(3, 2, 1));
          searching = false;
        }
        if (searching and compareedge(pi, pj, 3, 3)) {
          // std::cout << "   In some case3\n";
          out.push_back(setGlue(3, 3, -1));
          searching = false;
        }
      }
    }
  }

  out.shrink_to_fit();

  // std::cout << "I found " << out.size() << " patch interfaces.\n";
  // std::cout << "I will match\n";
  // for (auto g : out) {
  //   std::cout << "patch " << g.patch1 << " with Patch " << g.patch2 << "\n";
  //   std::cout << " edge " << g.matchp1 << " with edge  " << g.matchp2 <<
  //   "\n";
  // }

  return out;
}
}  // namespace Spl
#endif
