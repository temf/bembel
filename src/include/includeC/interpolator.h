// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_INTERPOLATION__
#define __BEMBEL_INTERPOLATION__

#include <Eigen/Sparse>
#include "Spline.hpp"
#include "glue.h"

namespace Bembel {
namespace Interpolation {
using _tpdata = std::array<int, 4>;

inline void stretch(Gluebucket &g, double const *in, double *out) {
  const int size = g.splinedofsnonglue;
  std::vector<int> tracker;
  tracker.reserve(g.size());

  auto isInTracker = [](int k, const std::vector<int> &v) {
    for (auto t : v) {
      if (t == k) return true;
    }
    return false;
  };
  int count = 0;
  // For every of the large dofs....
  for (int i = 0; i + count < size; i++) {
    // Checks if the dof has already been assigned a value
    while (isInTracker(i + count, tracker)) {
      count++;
    }
    // If not looks for a partner
    auto id = g.getPartnerId(i + count);
    // If no partner has been found, copy the next value of in
    if (id[0] == -1) {
      out[i + count] = in[i];
    } else {
      // If a partner has been found, copy the value to the next element of out
      // and the assigned parter up to sign
      out[i + count] = in[i];
      out[id[0]] = in[i] * id[1];
      // Track the element which has been assigned
      tracker.push_back(id[0]);
    }
  }
  return;
}

inline void stouch(Gluebucket &g, double const *in, double *out) {
  const int size = g.splinedofsnonglue;
  // The tracker tracks which large dofs have had an impact already
  std::vector<int> tracker;
  tracker.reserve(g.size());
  auto isInTracker = [](int k, const std::vector<int> &v) {
    for (auto t : v) {
      if (t == k) return true;
    }
    return false;
  };
  int count = 0;
  for (int i = 0; i + count < size; i++) {
    while (isInTracker(i + count, tracker)) {
      count++;
    }
    auto id = g.getPartnerId(i + count);
    // If no partner has been found, copy the next value of in
    if (id[0] == -1) {
      out[i] = in[i + count];
    } else {
      // If a partner has been found, copy the value to the next element of out
      // and the assigned parter up to sign
      out[i] = .5 * (in[i + count] + in[id[0]] * id[1]);
      // Track the element which has been assigned
      tracker.push_back(id[0]);
    }
  }
  return;
}

// This function figures out if the point pt is close to the support of the
// given dof within the spline space defined by p and knt.
inline bool isInSupport(const int dof, const int p,
                        const std::vector<double> &knt, const double pt) {
  constexpr double tol = 0.000000001;
  // const int dofs = knt.size() - p;
  return knt[dof] <= pt + tol and pt - tol <= knt[dof + p];
}

// dofisone evaluate dof 1 at the given point
inline double dofisOne(const int dof, const int deg,
                       const std::vector<double> &knt, const double pt) {
  // Eigen::Matrix<double,-1,1> coefx = Eigen::Matrix::Zeroes(kntx.size()-degx);
  Eigen::Matrix<double, -1, -1> coef =
      Eigen::MatrixXd::Zero(1, knt.size() - deg);
  // std::cout << "dofIsOne: 1," << knt.size() - deg << "\n";
  // std::cout << "dof is " << dof << "\n";
  // coefx(dofx) = 1;
  coef(dof) = 1;
  return BSPLN::deBoor(coef, knt, {pt})(0);
}

// This function takes in information about a uniformly refined TP-BSpline space
// (indat = degs,degt,refinement), as well as a pointer array (in) of
// coefficients. Same for out, where out* is a pointer array already of the
// correct size.
inline void interpolate(int length, _tpdata indat, double *in, _tpdata outdat,
                        double *out) {
  assert(indat.size() == 4 and outdat.size() == 4 and
         "Should be of size 4: (degs,degt,level,kntrep)");

  // constexpr double checkTol = 0.000001;

  // The number of interior element interface (interior knots up to repetition)
  const int outdegx = outdat[0];
  const int outdegy = outdat[1];
  const int indegx = indat[0];
  const int indegy = indat[1];
  const std::vector<double> outkntx =
      HLPR::mkUnifKnot(outdegx, outdat[2], outdat[3]);
  const std::vector<double> outknty =
      HLPR::mkUnifKnot(outdegy, outdat[2], outdat[3]);

  const std::vector<double> inkntx =
      HLPR::mkUnifKnot(indegx, indat[2], indat[3]);
  const std::vector<double> inknty =
      HLPR::mkUnifKnot(indegy, indat[2], indat[3]);
  const int in_s = inkntx.size() - indegx;
  const int in_t = inknty.size() - indegy;
  const int indofs = in_s * in_t;
  const int out_s = outkntx.size() - outdegx;
  const int out_t = outknty.size() - outdegy;
  const int outdofs = out_s * out_t;

  // std::cout << out_s << " and " << out_t << "\n";

  // std::cout << "Interpolator yields " << outdofs * length << " dofs.\n";

  // std::cout <<" outdosf = "<<outdofs << "\n";
  std::vector<std::array<double, 2>> intpts;
  const double factor = 1;
  // intpts.reserve(outdofs * 4);
  for (int ix = 0; ix < factor * out_s; ix++)
    for (int iy = 0; iy < factor * out_t; iy++)
      intpts.push_back({((double)(ix)) / (factor * out_s - 1),
                        ((double)(iy)) / (factor * out_t - 1)});

  for (auto x : intpts) {
    assert(-.00000000001 < x[0] and x[0] < 1.0000000001 and
           -.00000000001 < x[1] and x[1] < 1.0000000001);
  }

  // The coefficient matrices are dense; however, the inerpolation matrices are
  // not. I dont even need to recunstruc, since all patches are the same.

  // std::cout << "             Unit interpolation...     ";
  std::vector<Eigen::Triplet<double>> trips;
  // trips.reserve(outdofs * 25); // rough estimate....
  const int intptssize = intpts.size();

  {
    for (int intpt = 0; intpt < intptssize; intpt++) {
      for (int dofx = 0; dofx < out_s; dofx++) {
        if (isInSupport(dofx, outdegx, outkntx, intpts[intpt][0])) {
          for (int dofy = 0; dofy < out_t; dofy++) {
            if (isInSupport(dofy, outdegy, outknty, intpts[intpt][1])) {
              const double xval =
                  dofisOne(dofx, outdegx, outkntx, intpts[intpt][0]);
              const double yval =
                  dofisOne(dofy, outdegy, outknty, intpts[intpt][1]);
              trips.push_back(Eigen::Triplet<double>(intpt, dofx * out_t + dofy,
                                                     xval * yval));
            }
          }
        }
      }
    }
  }
  // std::cout << "...done.\n             Factorizing Sparse LU...  \n";

  Eigen::SparseMatrix<double, ColMajor> interpol(intpts.size(), outdofs);
  interpol.setFromTriplets(trips.begin(), trips.end());
  interpol.makeCompressed();

  SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
  solver.analyzePattern(interpol);
  solver.compute(interpol);

  // std::cout << "...done.\n             Evaluating the RHS...     \n";

  Eigen::Matrix<double, -1, 1> rhs(intpts.size());
  for (int k = 0; k < length; k++) {
    MatrixXd incoeffvec(in_t, in_s);

    for (int dofx = 0; dofx < in_s; dofx++)
      for (int dofy = 0; dofy < in_t; dofy++)
        incoeffvec(dofx * in_t + dofy) = in[k * indofs + dofx * in_t + dofy];

    // std::cout << incoeffvec << "\n";

    // for(auto x:inkntx)
    // std::cout << x << ",";
    // std::cout <<"\n";
    // for(auto x:inknty)
    // std::cout << x << ",";
    // std::cout <<"\n";
    for (int i = 0; i < intptssize; i++) {
      rhs(i) = BSPLN::deBoorTP(incoeffvec, inkntx, inknty, {intpts[i][0]},
                               {intpts[i][1]})(0);
    }

    // std::cout << "...done.\n             Solving...   ";

    Eigen::Matrix<double, -1, -1> x = solver.solve(rhs);
    for (int i = 0; i < outdofs; i++) out[k * outdofs + i] = x(i);

    // std::cout << "             ...done.\n";
    // assert((x(outdofs - 1) - incoeffvec(indofs - 1)) < checkTol &&
    //        (x(0) - incoeffvec(0)) < checkTol &&
    //        "This should be given due to the boundary interpolation property
    //        of " "any of the constructable spline spaces.");
  }

  // std::cout << "End interp.\n";
  // std::cout << "End interp.\n";
  // std::cout << "End interp.\n";
  // std::cout << "End interp.\n";
  return;
}

}  // namespace Interpolation
}  // namespace Bembel
#endif
