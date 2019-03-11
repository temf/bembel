// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "glue.h"

namespace Bembel {
// local helperstruct
struct Dofs {
  int xdir1comp;
  int ydir1comp;
  int xdir2comp;
  int ydir2comp;
  int patchdof1comp;
  int patchdof2comp;
  int componentdof1comp;
  int componentdof2comp;
};

Dofs dofmaker(int a, int b, int c, int d, int e, int f, int g, int h) {
  Dofs dof;
  dof.xdir1comp = a;
  dof.ydir1comp = b;
  dof.ydir2comp = c;
  dof.xdir2comp = d;
  dof.patchdof1comp = e;
  dof.patchdof2comp = f;
  dof.componentdof1comp = g;
  dof.componentdof2comp = h;
  return dof;
}

inline void printDofInfo(Dofs dof) {
  std::cout << "dof.xdir1comp = " << dof.xdir1comp << std::endl;
  std::cout << "dof.ydir1comp = " << dof.ydir1comp << std::endl;
  std::cout << "dof.ydir2comp = " << dof.ydir2comp << std::endl;
  std::cout << "dof.xdir2comp = " << dof.xdir2comp << std::endl;
  std::cout << "dof.patchdof1comp = " << dof.patchdof1comp << std::endl;
  std::cout << "dof.patchdof2comp = " << dof.patchdof2comp << std::endl;
  std::cout << "dof.componentdof1comp = " << dof.componentdof1comp << std::endl;
  std::cout << "dof.componentdof2comp = " << dof.componentdof2comp << std::endl;
}

// THIS IS FOR MAXWELL ONLY
// We need to glue ONCE per edge; and the numbers are generated w.r.t.
// reference patch numbering, i.e., starting at 0. This means that the
// output of this needs to be shifted to the correct place w.r.t. the
// projection.
// \brief This function gets informations about the edges that need to be glued
// together. Then, depending on the information given by glue g, it goes into
// switch statements which yield the IDs of Dofs that need to be glued together.
// It does not need the discretization for this: All it needs to know about the
// global space is that the numbering starts in y-direction first, and the
// numbers of 1D-dofs on which the TP constructino is based. All of this is
// initialized in the Dofs helper struct, which is (at the moment) only used
// here. I might make it part of the discretization settings-struct.
std::vector<std::array<int, 2>> getGlueDofIds_Maxwell(const Spl::glue &g,
                                                      const Dofs &dof) {
  // Now, we look at the dof numbers of the individual patch.
  // This switch does some arithmetic: The cases are :
  // 0 : x == 0 Kante
  // 1 : Y == 1 Kante
  // 2 : X == 1 Kante
  // 3 : Y == 0 Kante
  // See boundarymatch.h for more info.
  //
  // For now only with global a_o!

  const int smallPolyDeg = dof.ydir1comp;
  const int largePolyDeg = dof.xdir1comp;

  const int dofs_per_component_per_patch = dof.patchdof1comp;

  // Since (for now) the implementation should work with a uniform degree at the
  // start of the De Rham complex, i.e., degrees in different direction can only
  // differ by one, we check some numbers to check for mistakes.
  assert(dof.ydir1comp == dof.xdir2comp &&
         "These are the smaller polynomial "
         "degrees. They must match by "
         "assumption.");
  assert(dof.xdir1comp == dof.ydir2comp &&
         "These are the larger polynomial "
         "degrees. They must match by "
         "assumption.");
  assert((dof.xdir1comp - 1 == dof.ydir1comp) or
         (dof.xdir1comp + 1 == dof.ydir1comp) &&
             "The difference of polynomial degree must be one!");
  assert(dof.patchdof1comp == dof.patchdof2comp &&
         "Without identification of DOFS on interfaces, both vector components "
         "must have the same size.");
  assert(dof.patchdof1comp == dof.ydir1comp * dof.xdir1comp &&
         "The degrees of freedom per dof must be governed by a tensor product "
         "relation.");

  std::vector<int> fro, to;

  const int secndcompshift = dof.componentdof1comp;
  // Shifts to the dofs of second vector component
  const int complex1stcompshift = dof.componentdof1comp + dof.componentdof2comp;
  // Shifts to the complex dofs of first vector component
  const int complex2ndcompshift = complex1stcompshift + secndcompshift;
  // Shifts to the complex dofs of second vector component

  std::vector<std::array<int, 2>> out;
  ///////////////////////////////////////////////////////////

  switch (g.matchp1) {
    case (0):
      // Running along x starting from 0 are the smallPolyDeg first dofs
      for (int j = 0; j < smallPolyDeg; j++) {
        const int loc_dofid = j * largePolyDeg;
        const int loc_patchshift = g.patch1 * dofs_per_component_per_patch;
        fro.push_back(loc_dofid + loc_patchshift);
        fro.push_back(loc_dofid + loc_patchshift + complex1stcompshift);
      }
      break;
    // This is a Y edge, meaning the dofs are not right behind each other, i.e.,
    // arranged in increments of largePolyDeg. Moreover, we need to glue the
    // second vector component, and thus need to shift.
    case (1):
      for (int j = 0; j < smallPolyDeg; j++) {
        const int loc_dofid = smallPolyDeg * (largePolyDeg - 1) + j;
        const int loc_patchshift = g.patch1 * dofs_per_component_per_patch;
        fro.push_back(loc_dofid + loc_patchshift + secndcompshift);
        fro.push_back(loc_dofid + loc_patchshift + complex2ndcompshift);
      }
      break;
      // These are, in analogy to case 0, the last smallPolyDeg dofs on each
      // patch in the first vector component
    case (2):
      for (int j = 0; j < smallPolyDeg; j++) {
        const int loc_dofid = (j + 1) * largePolyDeg - 1;
        const int loc_patchshift = g.patch1 * dofs_per_component_per_patch;
        fro.push_back(loc_dofid + loc_patchshift);
        fro.push_back(loc_dofid + loc_patchshift + complex1stcompshift);
      }
      break;
    default:
      assert(g.matchp1 == 3);
      // This is essentially case 3. We need the last dof of each Dof-line in  X
      // direction. Moreover, we need to shift to the second vector component
      for (int j = 0; j < smallPolyDeg; j++) {
        const int loc_dofid = j;
        const int loc_patchshift = g.patch1 * dofs_per_component_per_patch;
        fro.push_back(loc_dofid + loc_patchshift + secndcompshift);
        fro.push_back(loc_dofid + loc_patchshift + complex2ndcompshift);
      }
      break;
  }

  ///////////////////////////////////////////////////////////
  switch (g.matchp2) {
    case (0):
      // Running along x starting from 0 are the smallPolyDeg first dofs
      for (int j = 0; j < smallPolyDeg; j++) {
        const int loc_dofid = j * largePolyDeg;
        const int loc_patchshift = g.patch2 * dofs_per_component_per_patch;
        if (g.reverse) {
          to.push_back(loc_dofid + loc_patchshift + complex1stcompshift);
        }
        to.push_back(loc_dofid + loc_patchshift);
        if (not g.reverse) {
          to.push_back(loc_dofid + loc_patchshift + complex1stcompshift);
        }
      }
      break;
    // This is a Y edge, meaning the dofs are not right behind each other, i.e.,
    // arranged in increments of largePolyDeg. Moreover, we need to glue the
    // second vector component, and thus need to shift.
    case (1):
      for (int j = 0; j < smallPolyDeg; j++) {
        const int loc_dofid = smallPolyDeg * (largePolyDeg - 1) + j;
        const int loc_patchshift = g.patch2 * dofs_per_component_per_patch;
        if (g.reverse) {
          to.push_back(loc_dofid + loc_patchshift + complex2ndcompshift);
        }
        to.push_back(loc_dofid + loc_patchshift + secndcompshift);
        if (not g.reverse) {
          to.push_back(loc_dofid + loc_patchshift + complex2ndcompshift);
        }
      }
      break;
      // These are, in analogy to case 0, the last smallPolyDeg dofs on each
      // patch in the first vector component
    case (2):
      for (int j = 0; j < smallPolyDeg; j++) {
        const int loc_dofid = (j + 1) * largePolyDeg - 1;
        const int loc_patchshift = g.patch2 * dofs_per_component_per_patch;
        if (g.reverse) {
          to.push_back(loc_dofid + loc_patchshift + complex1stcompshift);
        }
        to.push_back(loc_dofid + loc_patchshift);
        if (not g.reverse) {
          to.push_back(loc_dofid + loc_patchshift + complex1stcompshift);
        }
      }
      break;
    default:
      assert(g.matchp2 == 3);
      // This is essentially case 3. We need the last dof of each Dof-line in  X
      // direction. Moreover, we need to shift to the second vector component
      for (int j = 0; j < smallPolyDeg; j++) {
        const int loc_dofid = j;
        const int loc_patchshift = g.patch2 * dofs_per_component_per_patch;
        if (g.reverse) {
          to.push_back(loc_dofid + loc_patchshift + complex2ndcompshift);
        }
        to.push_back(loc_dofid + loc_patchshift + secndcompshift);
        if (not g.reverse) {
          to.push_back(loc_dofid + loc_patchshift + complex2ndcompshift);
        }
      }
      break;
  }

  const int sz = to.size();
  if (g.reverse) {
    std::reverse(to.begin(), to.end());
  }

  assert(to.size() == fro.size() && "That should match. glue.cpp");
  out.reserve(sz);
  for (int j = 0; j < sz; j++) {
    out.push_back({to[j], fro[j]});
  }
  return out;
}

/*
    \brief This function assembles a sparse Matrix that can be applied to the
   projector to force the normal continuity required by the divergence
   conforming space.
*/
Eigen::SparseMatrix<double> make_glue_matrix_maxwell(
    discretization *ds, int p, int kntrep, int M, const int maxsizeofmatrix) {
  auto gls = boundary_match(*ds->mesh->geom);
  // for (auto x:gls) {
  //   std::cout << x.patch1 << " " << x.patch2 << " " << x.matchp1 << " " <<
  //   x.matchp2 << " " << x.gluecond << std::endl;
  // }

  Gluebucket gluebucket(ds, gls, p, M, kntrep);

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
      count = count + 1;
    } else {
      trips.push_back(Eigen::Triplet<int>(i - count, i, 1));
      // This is negative if there is no partner!
      const std::array<int, 2> pid = gluebucket.getPartnerId(i);
      if (pid[0] >= 0) {
        takencareof.push_back(pid[0]);
        trips.push_back(Eigen::Triplet<int>(i - count, pid[0], pid[1]));
      }
    }
  }

  Eigen::SparseMatrix<double> gluemat(maxsizeofmatrix - count, maxsizeofmatrix);
  gluemat.setFromTriplets(trips.begin(), trips.end());

  return gluemat;
}

Gluebucket::Gluebucket(discretization *ds, const std::vector<Spl::glue> &glues,
                       int p, int M, int kntrep) {
  const int a_o = ds->a_o;
  const int largeX = (a_o) + ((1 << M) - 1) * std::min(a_o, kntrep);
  const int largeY = (a_o) + ((1 << M) - 1) * std::min(a_o, kntrep);
  const int smallX = (a_o - 1) + ((1 << M) - 1) * std::min(a_o - 1, kntrep);
  const int smallY = (a_o - 1) + ((1 << M) - 1) * std::min(a_o - 1, kntrep);
  const int patchn = ds->mesh->geom->size();

  const int comp1 = largeX * smallX;
  const int comp2 = largeY * smallY;
  Dofs dofs = dofmaker(largeX, smallX, largeY, smallY, comp1, comp2,
                       patchn * comp1, patchn * comp2);

  // I assume that the order of the basis functions is:
  // realX,complexX,realY,complexY

  for (auto g : glues) {
    std::vector<std::array<int, 2>> ids = getGlueDofIds_Maxwell(g, dofs);
    for (auto t : ids) {
      if (not(t[0] == t[1])) {
        // gluefrom.push_back(std::min(t[0], t[1]));
        // glueto.push_back(std::max(t[0], t[1]));
        gluefrom.push_back(t[1]);
        glueto.push_back(t[0]);
        directions.push_back(g.gluecond);
      }
    }
  }

  splinedofsnonglue = 2 * patchn * (comp1 + comp2);

  gluefrom.shrink_to_fit();
  glueto.shrink_to_fit();
  directions.shrink_to_fit();
}
}  // namespace Bembel