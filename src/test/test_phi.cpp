// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "bemtest.h"

using namespace Bembel;

int Test::test_phi() {
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

  for (int iorder = 2; iorder < 11; ++iorder) {
    int iorderp1 = iorder + 1;

    // get coefficients for interpolation of f and g
    double t[iorder];
    for (int i = 0; i < iorder; ++i) t[i] = (double)i / (iorder - 1);

    // Here, we suddenly use the degree for basisevaluation, i.e., maxp-1. This
    // is confusing, but correct and tested.
    Eigen::Matrix<double, -1, -1> system(iorder, iorder);
    for (int iy = 0; iy < iorder; iy++) {
      std::vector<double> vals_y = Spl::evalBrnstnBasis(iorder - 1, t[iy]);
      for (int jy = 0; jy < iorder; jy++) {
        system(iy, jy) = vals_y[jy];
      }
    }

    // Do LU-Decomposition of system
    const Eigen::PartialPivLU<Eigen::Matrix<double, -1, -1>> pplu(system);
    // assemble and solve interpolation problems
    Eigen::VectorXd rhsf(iorder);
    Eigen::VectorXd rhsg(iorder);
    for (int iy = 0; iy < iorder; iy++) {
      rhsf(iy) = f(t[iy], iorder);
      rhsg(iy) = g(t[iy], iorder);
    }
    Eigen::VectorXd solf = pplu.solve(rhsf);
    Eigen::VectorXd solg = pplu.solve(rhsg);

    // get rid of Eigen
    double f_coeff[iorder];
    double g_coeff[iorder];
    for (int i = 0; i < iorder; ++i) {
      f_coeff[i] = solf(i);
      g_coeff[i] = solg(i);
    }

    // check whether coefficients are correct
    {
      double t[2 * iorder];
      for (int i = 0; i < 2 * iorder; ++i) t[i] = (double)i / (2 * iorder - 1);
      double testf[iorder];
      double testg[iorder];
      memset(testf, 0, iorder * sizeof(double));
      memset(testg, 0, iorder * sizeof(double));
      for (int i = 0; i < iorder; ++i) {
        double test[iorder];
        memset(test, 0, iorder * sizeof(double));
        select_phi(iorder)(test, 1., t[i]);
        testf[i] = myddot(iorder, f_coeff, test);
        testg[i] = myddot(iorder, g_coeff, test);
      }
      for (int i = 0; i < iorder; ++i) {
        bool_assert(fabs(testf[i] - f(t[i], iorder)) < 1e-14);
        bool_assert(fabs(testg[i] - g(t[i], iorder)) < 1e-14);
      }
    }

    // check phi_dx
    {
      double t[2 * iorder];
      for (int i = 0; i < 2 * iorder; ++i) t[i] = (double)i / (2 * iorder - 1);
      double testf[iorder];
      double testg[iorder];
      memset(testf, 0, iorder * sizeof(double));
      memset(testg, 0, iorder * sizeof(double));
      for (int i = 0; i < iorder; ++i) {
        double test[iorder];
        memset(test, 0, iorder * sizeof(double));
        select_phi_dx(iorder)(test, 1., t[i]);
        testf[i] = myddot(iorder, f_coeff, test);
        testg[i] = myddot(iorder, g_coeff, test);
      }
      for (int i = 0; i < iorder; ++i) {
        bool_assert(fabs(testf[i] - f_dx(t[i], iorder)) < 1e-14);
        bool_assert(fabs(testg[i] - g_dx(t[i], iorder)) < 1e-14);
      }
    }

    // std::cout << "Order " << iorder << " is ok." << std::endl;
  }

  return (0);
}
