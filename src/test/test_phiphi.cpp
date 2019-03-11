// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "bemtest.h"

using namespace Bembel;
int Test::test_phiphi() {
  int count = 0;

  auto bool_assert = [&count](bool in) {
    if (in) {
      return;
    } else
      count++;
    return;
  };

  auto f = [](double x, int iorder) { return pow(x, iorder - 1); };
  auto f_dx = [](double x, int iorder) {
    return (iorder - 1) * pow(x, iorder - 2);
  };
  auto g = [](double x, int iorder) { return pow(1. - x, iorder - 1); };
  auto g_dx = [](double x, int iorder) {
    return -(iorder - 1) * pow(1. - x, iorder - 2);
  };

  auto fg = [&f, &g](vector2 a, int iorder) {
    return f(a.x, iorder) * g(a.y, iorder);
  };
  auto fg_dx = [&f_dx, &g](vector2 a, int iorder) {
    return f_dx(a.x, iorder) * g(a.y, iorder);
  };
  auto fg_dy = [&f, &g_dx](vector2 a, int iorder) {
    return f(a.x, iorder) * g_dx(a.y, iorder);
  };

  for (int iorder = 2; iorder < 11; ++iorder) {
    int iorder2 = iorder * iorder;

    // get coefficients for interpolation of f and g
    vector2 t[iorder2];
    for (int i = 0; i < iorder; ++i)
      for (int j = 0; j < iorder; ++j)
        t[i * iorder + j] =
            vector2_make((double)j / (iorder - 1), (double)i / (iorder - 1));

    // Here, we suddenly use the degree for basisevaluation, i.e., maxp-1. This
    // is confusing, but correct and tested.
    Eigen::Matrix<double, -1, -1> system(iorder2, iorder2);
    for (int iy = 0; iy < iorder; iy++) {
      std::vector<double> vals_y =
          Spl::evalBrnstnBasis(iorder - 1, t[iy * iorder].y);
      for (int ix = 0; ix < iorder; ix++) {
        std::vector<double> vals_x =
            Spl::evalBrnstnBasis(iorder - 1, t[iy * iorder + ix].x);
        for (int jy = 0; jy < iorder; jy++) {
          for (int jx = 0; jx < iorder; jx++) {
            system(iy * iorder + ix, jy * iorder + jx) =
                vals_x[jx] * vals_y[jy];
          }
        }
      }
    }
    // Do LU-Decomposition of system
    const Eigen::PartialPivLU<Eigen::Matrix<double, -1, -1>> pplu(system);
    // assemble and solve interpolation problems
    Eigen::VectorXd rhs(iorder2);
    for (int iy = 0; iy < iorder; iy++)
      for (int ix = 0; ix < iorder; ix++)
        rhs(iy * iorder + ix) = fg(t[iy * iorder + ix], iorder);
    Eigen::VectorXd sol = pplu.solve(rhs);

    // get rid of Eigen
    double fg_coeff[iorder2];
    for (int i = 0; i < iorder2; ++i) fg_coeff[i] = sol(i);

    // check whether coefficients are correct
    {
      vector2 t[4 * iorder2];
      for (int i = 0; i < iorder; ++i)
        for (int j = 0; j < iorder; ++j)
          t[i * iorder + j] = vector2_make((double)j / (2 * iorder - 1),
                                           (double)i / (2 * iorder - 1));
      double testf[iorder2];
      memset(testf, 0, iorder2 * sizeof(double));
      for (int i = 0; i < iorder2; ++i) {
        double test[iorder2];
        memset(test, 0, iorder2 * sizeof(double));
        select_phiphi(iorder)(test, t[i]);
        testf[i] = myddot(iorder2, fg_coeff, test);
      }
      for (int i = 0; i < iorder2; ++i) {
        bool_assert(fabs(testf[i] - fg(t[i], iorder)) < 1e-14);
      }
    }

    // check whether coefficients are correct
    {
      vector2 t[4 * iorder2];
      for (int i = 0; i < iorder; ++i)
        for (int j = 0; j < iorder; ++j)
          t[i * iorder + j] = vector2_make((double)j / (2 * iorder - 1),
                                           (double)i / (2 * iorder - 1));
      double testf[iorder2];
      memset(testf, 0, iorder2 * sizeof(double));
      for (int i = 0; i < iorder2; ++i) {
        double test[iorder2];
        memset(test, 0, iorder2 * sizeof(double));
        select_phiphi_dx(iorder)(test, t[i]);
        testf[i] = myddot(iorder2, fg_coeff, test);
      }
      for (int i = 0; i < iorder2; ++i) {
        bool_assert(fabs(testf[i] - fg_dx(t[i], iorder)) < 1e-14);
      }
    }

    // check whether coefficients are correct
    {
      vector2 t[4 * iorder2];
      for (int i = 0; i < iorder; ++i)
        for (int j = 0; j < iorder; ++j)
          t[i * iorder + j] = vector2_make((double)j / (2 * iorder - 1),
                                           (double)i / (2 * iorder - 1));
      double testf[iorder2];
      memset(testf, 0, iorder2 * sizeof(double));
      for (int i = 0; i < iorder2; ++i) {
        double test[iorder2];
        memset(test, 0, iorder2 * sizeof(double));
        select_phiphi_dy(iorder)(test, t[i]);
        testf[i] = myddot(iorder2, fg_coeff, test);
      }
      for (int i = 0; i < iorder2; ++i) {
        bool_assert(fabs(testf[i] - fg_dy(t[i], iorder)) < 1e-14);
      }
    }

    // std::cout << "Order " << iorder << " is ok." << std::endl;
  }

  return count;
}
