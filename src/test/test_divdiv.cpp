// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "bemtest.h"

using namespace Bembel;

int Test::test_divdiv() {
  int count = 0;

  auto bool_assert = [&count](bool in) {
    if (in) {
      return;
    } else
      count++;
    return;
  };

  auto f_1 = [](vector2 a) { return a.x * a.y; };
  auto f_2 = [](vector2 a) { return -a.x * a.y; };
  auto f_1_dx = [](vector2 a) { return a.y; };
  auto f_2_dy = [](vector2 a) { return -a.x; };

  auto g_1 = [](vector2 a) { return -a.x * a.y; };
  auto g_2 = [](vector2 a) { return a.x * a.y; };
  auto g_1_dx = [](vector2 a) { return -a.y; };
  auto g_2_dy = [](vector2 a) { return a.x; };

  auto div_f1g1 = [&f_1_dx, &g_1_dx](vector2 a, vector2 b) {
    return f_1_dx(a) * g_1_dx(b);
  };
  auto div_f1g2 = [&f_1_dx, &g_2_dy](vector2 a, vector2 b) {
    return f_1_dx(a) * g_2_dy(b);
  };
  auto div_f2g1 = [&f_2_dy, &g_1_dx](vector2 a, vector2 b) {
    return f_2_dy(a) * g_1_dx(b);
  };
  auto div_f2g2 = [&f_2_dy, &g_2_dy](vector2 a, vector2 b) {
    return f_2_dy(a) * g_2_dy(b);
  };

  for (int iorder = 2; iorder < 11; ++iorder) {
    int iorderp1 = iorder + 1;
    int iorder2 = iorder * iorder;
    int iorderp12 = iorderp1 * iorderp1;

    // get coefficients for interpolation of f and g
    double f1[iorder2];
    double f2[iorder2];
    double g1[iorder2];
    double g2[iorder2];
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
    Eigen::VectorXd rhsf1(iorder2);
    Eigen::VectorXd rhsf2(iorder2);
    Eigen::VectorXd rhsg1(iorder2);
    Eigen::VectorXd rhsg2(iorder2);
    for (int iy = 0; iy < iorder; iy++)
      for (int ix = 0; ix < iorder; ix++) {
        rhsf1(iy * iorder + ix) = f_1(t[iy * iorder + ix]);
        rhsf2(iy * iorder + ix) = f_2(t[iy * iorder + ix]);
        rhsg1(iy * iorder + ix) = g_1(t[iy * iorder + ix]);
        rhsg2(iy * iorder + ix) = g_2(t[iy * iorder + ix]);
      }
    Eigen::VectorXd solf1 = pplu.solve(rhsf1);
    Eigen::VectorXd solf2 = pplu.solve(rhsf2);
    Eigen::VectorXd solg1 = pplu.solve(rhsg1);
    Eigen::VectorXd solg2 = pplu.solve(rhsg2);
    for (int i = 0; i < iorder; ++i)
      for (int j = 0; j < iorder; ++j) {
        f1[i * iorder + j] = solf1(i * iorder + j);
        f2[i * iorder + j] = solf2(i * iorder + j);
        g1[i * iorder + j] = solg1(i * iorder + j);
        g2[i * iorder + j] = solg2(i * iorder + j);
      }

    // new evaluation points
    for (int i = 0; i < iorder; ++i)
      for (int j = 0; j < iorder; ++j)
        t[i * iorder + j] = vector2_make((double)(j + 1) / (iorder + 1),
                                         (double)(i + 1) / (iorder + 1));

    // check whether coefficients are correct
    double testf1[iorder2];
    double testf2[iorder2];
    double testg1[iorder2];
    double testg2[iorder2];
    memset(testf1, 0, iorder2 * sizeof(double));
    memset(testf2, 0, iorder2 * sizeof(double));
    memset(testg1, 0, iorder2 * sizeof(double));
    memset(testg2, 0, iorder2 * sizeof(double));
    {
      for (int i = 0; i < iorder2; ++i) {
        double test[iorder2];
        memset(test, 0, iorder2 * sizeof(double));
        select_phiphi(iorder)(test, t[i]);
        testf1[i] = myddot(iorder2, f1, test);
        testf2[i] = myddot(iorder2, f2, test);
        testg1[i] = myddot(iorder2, g1, test);
        testg2[i] = myddot(iorder2, g2, test);
      }
      for (int i = 0; i < iorder2; ++i) {
        bool_assert(fabs(testf1[i] - f_1(t[i])) < 1e-12);
        bool_assert(fabs(testf2[i] - f_2(t[i])) < 1e-12);
        bool_assert(fabs(testg1[i] - g_1(t[i])) < 1e-12);
        bool_assert(fabs(testg2[i] - g_2(t[i])) < 1e-12);
      }
    }

    // check whether coefficients are also correct for the derivatives
    double testf1_dx[iorder2];
    double testf2_dy[iorder2];
    double testg1_dx[iorder2];
    double testg2_dy[iorder2];
    memset(testf1_dx, 0, iorder2 * sizeof(double));
    memset(testf2_dy, 0, iorder2 * sizeof(double));
    memset(testg1_dx, 0, iorder2 * sizeof(double));
    memset(testg2_dy, 0, iorder2 * sizeof(double));
    {
      for (int i = 0; i < iorder2; ++i) {
        double test[iorder2];
        memset(test, 0, iorder2 * sizeof(double));
        select_phiphi_dx(iorder)(test, t[i]);
        testf1_dx[i] = myddot(iorder2, f1, test);
        testg1_dx[i] = myddot(iorder2, g1, test);
        memset(test, 0, iorder2 * sizeof(double));
        select_phiphi_dy(iorder)(test, t[i]);
        testf2_dy[i] = myddot(iorder2, f2, test);
        testg2_dy[i] = myddot(iorder2, g2, test);
      }
      for (int i = 0; i < iorder2; ++i) {
        bool_assert(fabs(testf1_dx[i] - f_1_dx(t[i])) < 1e-12);
        bool_assert(fabs(testf2_dy[i] - f_2_dy(t[i])) < 1e-12);
        bool_assert(fabs(testg1_dx[i] - g_1_dx(t[i])) < 1e-12);
        bool_assert(fabs(testg2_dy[i] - g_2_dy(t[i])) < 1e-12);
      }
    }

    // check whether VPhi_scal_Phi is correct
    {
      double cf1g1[iorder2 * iorder2];
      double cf1g2[iorder2 * iorder2];
      double cf2g1[iorder2 * iorder2];
      double cf2g2[iorder2 * iorder2];
      for (int i = 0; i < iorder2; ++i)
        for (int j = 0; j < iorder2; ++j) {
          double c[4 * iorder2 * iorder2];
          double ct1[iorder2];
          double ct2[iorder2];
          double ct3[iorder2];
          double ct4[iorder2];
          memset(c, 0, 4 * iorder2 * iorder2 * sizeof(double));
          select_VPhi_scal_bla(iorder)(
              c, 1., t[i], t[j], vector3_make(1., 0., 0.),
              vector3_make(0., 1., 0.), vector3_make(1., 0., 0.),
              vector3_make(0., 1., 0.));
          for (int k = 0; k < iorder2; ++k) {
            ct1[k] =
                myddot(iorder2, c + 0 * iorder2 * iorder2 + k * iorder2, g1);
            ct2[k] =
                myddot(iorder2, c + 1 * iorder2 * iorder2 + k * iorder2, g2);
            ct3[k] =
                myddot(iorder2, c + 2 * iorder2 * iorder2 + k * iorder2, g1);
            ct4[k] =
                myddot(iorder2, c + 3 * iorder2 * iorder2 + k * iorder2, g2);
          }
          cf1g1[i * iorder2 + j] = myddot(iorder2, ct1, f1);
          cf1g2[i * iorder2 + j] = myddot(iorder2, ct2, f1);
          cf2g1[i * iorder2 + j] = myddot(iorder2, ct3, f2);
          cf2g2[i * iorder2 + j] = myddot(iorder2, ct4, f2);
        }

#if 0
      // std::cout << "Reference solution is:" << std::endl;
      printf("%8.6g | ", 0.);
      for (int j = 0 ; j < iorder2 ; ++j)
        printf("%8.6g ", g_1(t[j]));
      for (int j = 0 ; j < iorder2 ; ++j)
        printf("%8.6g ", g_2(t[j]));
      // std::cout << std::endl;
      // std::cout << "----------------------------------------" << std::endl;
      for (int i = 0 ; i < iorder2 ; ++i) {
        printf("%8.6g | ", f_1(t[i]));
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", f_1(t[i])*g_1(t[j]));
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", 0.);
        // std::cout << std::endl;
      }
      for (int i = 0 ; i < iorder2 ; ++i) {
        printf("%8.6g | ", f_2(t[i]));
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", 0.);
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", f_2(t[i])*g_2(t[j]));
        // std::cout << std::endl;
      }
      // std::cout << std::endl;
      // std::cout << "My solution is:" << std::endl;
      printf("%8.6g | ", 0.);
      for (int j = 0 ; j < iorder2 ; ++j)
        printf("%8.6g ", testg1[j]);
      for (int j = 0 ; j < iorder2 ; ++j)
        printf("%8.6g ", testg2[j]);
      // std::cout << std::endl;
      // std::cout << "----------------------------------------" << std::endl;
      for (int i = 0 ; i < iorder2 ; ++i) {
        printf("%8.6g | ", testf1[i]);
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", cf1g1[i*iorder2+j]);
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", cf1g2[i*iorder2+j]);
        // std::cout << std::endl;
      }
      for (int i = 0 ; i < iorder2 ; ++i) {
        printf("%8.6g | ", testf2[i]);
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", cf2g1[i*iorder2+j]);
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", cf2g2[i*iorder2+j]);
        // std::cout << std::endl;
      }
      // std::cout << std::endl;
#endif

      for (int i = 0; i < iorder2; ++i)
        for (int j = 0; j < iorder2; ++j) {
          bool_assert(fabs(cf1g1[i * iorder2 + j] - f_1(t[i]) * g_1(t[j])) <
                      1e-12);
          bool_assert(fabs(cf1g2[i * iorder2 + j]) < 1e-12);
          bool_assert(fabs(cf2g1[i * iorder2 + j]) < 1e-12);
          bool_assert(fabs(cf2g2[i * iorder2 + j] - f_2(t[i]) * g_2(t[j])) <
                      1e-12);
        }
    }

    // check whether Div_Phi_scal_Div_Phi is correct
    // this is the same as above, but slightly modified
    {
      double cf1g1[iorder2 * iorder2];
      double cf1g2[iorder2 * iorder2];
      double cf2g1[iorder2 * iorder2];
      double cf2g2[iorder2 * iorder2];
      for (int i = 0; i < iorder2; ++i)
        for (int j = 0; j < iorder2; ++j) {
          double c[4 * iorder2 * iorder2];
          double ct1[iorder2];
          double ct2[iorder2];
          double ct3[iorder2];
          double ct4[iorder2];
          memset(c, 0, 4 * iorder2 * iorder2 * sizeof(double));
          select_Div_Phi_times_Div_Phi(iorder)(c, 1., t[i], t[j]);
          for (int k = 0; k < iorder2; ++k) {
            ct1[k] =
                myddot(iorder2, c + 0 * iorder2 * iorder2 + k * iorder2, g1);
            ct2[k] =
                myddot(iorder2, c + 1 * iorder2 * iorder2 + k * iorder2, g2);
            ct3[k] =
                myddot(iorder2, c + 2 * iorder2 * iorder2 + k * iorder2, g1);
            ct4[k] =
                myddot(iorder2, c + 3 * iorder2 * iorder2 + k * iorder2, g2);
          }
          cf1g1[i * iorder2 + j] = myddot(iorder2, ct1, f1);
          cf1g2[i * iorder2 + j] = myddot(iorder2, ct2, f1);
          cf2g1[i * iorder2 + j] = myddot(iorder2, ct3, f2);
          cf2g2[i * iorder2 + j] = myddot(iorder2, ct4, f2);
        }

#if 0
      // std::cout << "Reference solution is:" << std::endl;
      printf("%8.6g | ", 0.);
      for (int j = 0 ; j < iorder2 ; ++j)
        printf("%8.6g ", g_1_dx(t[j]));
      for (int j = 0 ; j < iorder2 ; ++j)
        printf("%8.6g ", g_2_dy(t[j]));
      // std::cout << std::endl;
      // std::cout << "----------------------------------------" << std::endl;
      for (int i = 0 ; i < iorder2 ; ++i) {
        printf("%8.6g | ", f_1_dx(t[i]));
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", div_f1g1(t[i],t[j]));
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", div_f1g2(t[i],t[j]));
        // std::cout << std::endl;
      }
      for (int i = 0 ; i < iorder2 ; ++i) {
        printf("%8.6g | ", f_2_dy(t[i]));
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", div_f2g1(t[i],t[j]));
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", div_f2g2(t[i],t[j]));
        // std::cout << std::endl;
      }
      // std::cout << std::endl;
      // std::cout << "My solution is:" << std::endl;
      printf("%8.6g | ", 0.);
      for (int j = 0 ; j < iorder2 ; ++j)
        printf("%8.6g ", testg1_dx[j]);
      for (int j = 0 ; j < iorder2 ; ++j)
        printf("%8.6g ", testg2_dy[j]);
      // std::cout << std::endl;
      // std::cout << "----------------------------------------" << std::endl;
      for (int i = 0 ; i < iorder2 ; ++i) {
        printf("%8.6g | ", testf1_dx[i]);
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", cf1g1[i*iorder2+j]);
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", cf1g2[i*iorder2+j]);
        // std::cout << std::endl;
      }
      for (int i = 0 ; i < iorder2 ; ++i) {
        printf("%8.6g | ", testf2_dy[i]);
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", cf2g1[i*iorder2+j]);
        for (int j = 0 ; j < iorder2 ; ++j)
          printf("%8.6g ", cf2g2[i*iorder2+j]);
        // std::cout << std::endl;
      }
      // std::cout << std::endl;
#endif

      for (int i = 0; i < iorder2; ++i)
        for (int j = 0; j < iorder2; ++j) {
          bool_assert(fabs(cf1g1[i * iorder2 + j] - div_f1g1(t[i], t[j])) <
                      1e-12);
          bool_assert(fabs(cf1g2[i * iorder2 + j] - div_f1g2(t[i], t[j])) <
                      1e-12);
          bool_assert(fabs(cf2g1[i * iorder2 + j] - div_f2g1(t[i], t[j])) <
                      1e-12);
          bool_assert(fabs(cf2g2[i * iorder2 + j] - div_f2g2(t[i], t[j])) <
                      1e-12);
        }
    }

    // std::cout << "Order " << iorder << " is ok." << std::endl;
  }

  return count;
}
