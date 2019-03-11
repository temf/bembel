// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "bemtest.h"

using namespace Bembel;

/**
 *  \brief         This functions tests the projector against the
 * Cox-DeBoor-formula
 *
 */
int Test::test_projector() {
  int out = 0;

  // Change this to add additional cases
  std::vector<int> orders = {2, 3};

  for (int ii = 0; ii < 2; ++ii) {
    int iorder = orders[ii];
    // pde = pdeproblem_LaplaceSingle();
    // pde.kappa[0] = 0;
    // pde.kappa[1] = 0;

    auto f = [&](double x, double y) {
// 1 -> Spline and 0->polynomial
#if 0 
// 1 -> knotrepetition 0-> not
#if 0
      std::vector<double> testknt = {0,0,0,.5,.5,1,1,1};
      Eigen::MatrixXd ctrl(5,5);
      ctrl <<  0,1,2,-2,3,//
               -1,1,1,5,-1,//
               -1,1,1,-5,-1,//
               -1,1,1,-5,-1,//                              
               3,2,1,0,0;//
#else
      std::vector<double> testknt = {0, 0, 0, .5, 1, 1, 1};
      Eigen::MatrixXd ctrl(4, 4);
      ctrl << 0, 1, 2, -2, //
          -1, 1, 1, -1,    //
          -1, 1, 1, -5,    //
          3, 2, 1, 0, 0;   //
#endif
      returSpldeBoorTP(ctrl, testknt, testknt, {x}, {y})(0);
#else
      return std::pow(x, iorder) * std::pow(1. - y, iorder);
#endif
    };

    for (int M = 2; M < 4; M++) {
      const int n = (1 << M);
      const int kntrep = 2;
      /*
       * initialize geometry
       */
      Geometry geom(Bembel::Test::mkScreen());

      /*
       * initialize discretization
       */

      LaplaceSingle laplace;

      Discretization<LaplaceSingle> modern_disc;
      modern_disc.init_Discretization(geom, laplace, iorder, kntrep, M);
      discretization disc = modern_disc.get_disc();

      std::vector<double> xs;

      std::vector<double> mask = Spl::make_interpolation_mask(iorder + 1);

      et_node *E = disc.mesh->E.patch[0];

      // find lowest level in element tree
      while (E->son[0]) E = E->son[0];

      const int syssize = n * n * (iorder + 1) * (iorder + 1);
      Eigen::VectorXd rhs_big(syssize);
      // assembling the rhs in row-major w.r.t. x

      std::vector<std::array<double, 2>> phys_dom_big_pt;
      std::vector<std::array<double, 2>> ref_dom_big_pt;

      for (int element = 0; element < n * n; element++) {
        for (int iy = 0; iy < iorder + 1; iy++) {
          for (int ix = 0; ix < iorder + 1; ix++) {
            const double pos_x = (double)E[element].index_s / n;
            const double pos_y = (double)E[element].index_t / n;
            phys_dom_big_pt.push_back(
                {pos_x + (1. / n * mask[ix]), pos_y + (1. / n * mask[iy])});
            ref_dom_big_pt.push_back({mask[ix], mask[iy]});
            rhs_big((element * (iorder + 1) * (iorder + 1)) +
                    (iy * (iorder + 1)) + ix) =
                f(pos_x + (1. / n * mask[ix]), pos_y + (1. / n * mask[iy]));
          }
        }
      }

      Eigen::Matrix<double, -1, -1> system_big =
          Eigen::Matrix<double, -1, -1>::Zero(syssize, syssize);

      const int phys_dom_big_pt_sz = phys_dom_big_pt.size();
      for (int pt = 0; pt < phys_dom_big_pt_sz; pt++) {
        double *pp_vals;
        pp_vals = (double *)calloc(sizeof(double), (iorder + 1) * (iorder + 1));
        disc.phiphi(pp_vals,
                    vector2_make(ref_dom_big_pt[pt][0], ref_dom_big_pt[pt][1]));
        for (int el = 0; el < n * n; el++) {
          const double x_min = (double)E[el].index_s / n;
          const double x_max = (double)E[el].index_s / n + (1. / n);
          const double y_min = (double)E[el].index_t / n;
          const double y_max = (double)E[el].index_t / n + (1. / n);
          const double x = phys_dom_big_pt[pt][0];
          const double y = phys_dom_big_pt[pt][1];
          if (x_min <= x and x <= x_max and y_min <= y and y <= y_max) {
            for (int el_dof = 0; el_dof < (iorder + 1) * (iorder + 1);
                 el_dof++) {
              system_big(pt, el * (iorder + 1) * (iorder + 1) + el_dof) =
                  pp_vals[el_dof];
            }
          }
        }
        free(pp_vals);
      }

      const Eigen::FullPivLU<Eigen::Matrix<double, -1, -1>> pplu_big(
          system_big);

      Eigen::Matrix<double, -1, -1> sol_big = pplu_big.solve(rhs_big);

      // Setting up the projector

      const auto info =
          init_projector_base(&disc, iorder + 1, iorder + 1, kntrep);

      const int size = info.vals.size();
      std::vector<Eigen::Triplet<double>> trips;

      for (int k = 0; k < size; k++) {
        trips.push_back(
            Eigen::Triplet<double>(info.rows[k], info.cols[k], info.vals[k]));
      }

      Eigen::SparseMatrix<double> Proj(info.num_dc_dof, info.num_c_dof);
      Proj.setFromTriplets(trips.begin(), trips.end());

      Eigen::MatrixXd proj = Proj;
      // std::cout << proj << std::endl << "\n";

      // knots etc.... for the continuous space
      std::vector<double> knot_x =
          Spl::make_unif_knots(iorder + 1, n - 1, kntrep);
      std::vector<double> knot_y =
          Spl::make_unif_knots(iorder + 1, n - 1, kntrep);
      const int num_c_dof_y = knot_y.size() - (iorder + 1);
      const int num_c_dof_x = knot_x.size() - (iorder + 1);
#ifdef _BEMBEL_TEST_WITH_ASSERT_
      assert(info.num_c_dof == num_c_dof_x * num_c_dof_y);
      assert(info.num_dc_dof == n * n * (iorder + 1) * (iorder + 1));
#endif
      if (not(info.num_c_dof == num_c_dof_x * num_c_dof_y) and
          (info.num_dc_dof == n * n * (iorder + 1) * (iorder + 1))) {
        out = 1;
      }

      Eigen::VectorXd rhs_small(num_c_dof_y * num_c_dof_x);
      std::vector<std::array<double, 2>> ref_dom_sma_pt;

      {
        std::vector<double> xs_s(num_c_dof_x);
        std::vector<double> ys_s(num_c_dof_y);

        for (int i = 0; i < num_c_dof_x; i++)
          xs_s[i] = i * (1. / (num_c_dof_x - 1));
        for (int i = 0; i < num_c_dof_y; i++)
          ys_s[i] = i * (1. / (num_c_dof_y - 1));
        for (int c_dof_y = 0; c_dof_y < num_c_dof_y; c_dof_y++) {
          for (int c_dof_x = 0; c_dof_x < num_c_dof_x; c_dof_x++) {
            ref_dom_sma_pt.push_back({xs_s[c_dof_x], ys_s[c_dof_y]});
            rhs_small(c_dof_x + c_dof_y * num_c_dof_x) =
                f(xs_s[c_dof_x], ys_s[c_dof_y]);
          }
        }
      }

      Eigen::Matrix<double, -1, -1> system_small(info.num_c_dof,
                                                 info.num_c_dof);

      {
        const int ref_dom_sma_pt_sz = ref_dom_sma_pt.size();
        std::vector<double> coef_x(num_c_dof_x, 0.0);
        std::vector<double> coef_y(num_c_dof_y, 0.0);
        for (int c_dof_y = 0; c_dof_y < num_c_dof_y; c_dof_y++) {
          coef_y[c_dof_y] = 1;
          for (int c_dof_x = 0; c_dof_x < num_c_dof_x; c_dof_x++) {
            coef_x[c_dof_x] = 1;
            for (int i = 0; i < ref_dom_sma_pt_sz; i++) {
              double x = Spl::deBoor(coef_x, knot_x, {ref_dom_sma_pt[i][0]})[0];
              double y = Spl::deBoor(coef_y, knot_y, {ref_dom_sma_pt[i][1]})[0];
              system_small(i, c_dof_y * num_c_dof_x + c_dof_x) = x * y;
            }
            coef_x[c_dof_x] = 0;
          }
          coef_y[c_dof_y] = 0;
        }
      }
      // std::cout << "big\n" << system_big << "\n" << std::endl;

      // std::cout << "small\n" << system_small << "\n" << std::endl;

      const Eigen::FullPivLU<Eigen::Matrix<double, -1, -1>> pplu_small(
          system_small);
#ifdef _BEMBEL_TEST_WITH_ASSERT_
      assert(pplu_small.rank() == info.num_c_dof);
      assert(pplu_big.rank() == info.num_dc_dof);
#endif
      Eigen::Matrix<double, -1, -1> sol_small = pplu_small.solve(rhs_small);

      Eigen::Matrix<double, -1, -1> err = Proj * sol_small - sol_big;

#ifdef _BEMBEL_PRINT_INFO_
      std::cout << "                  iorder : " << iorder << std::endl;
      std::cout << "                       M : " << M << std::endl;
      std::cout << "          Projector-Error: " << err.norm() << std::endl;
#endif

      if (err.norm() > 0.000001) {
        out = 1;
        Eigen::Matrix<double, -1, 3> tmp(err.rows(), 3);
        tmp.col(0) = (Proj * sol_small);
        tmp.col(1) = sol_big;
        tmp.col(2) = err;
        std::cout << "\n" << tmp << std::endl;
        // std::cout << (Proj * sol_small).transpose() << std::endl
        //           << "\n"
        //           << sol_big.transpose() << std::endl;
      }
#ifdef _BEMBEL_TEST_WITH_ASSERT_
      assert(err.norm() < 0.00001);
#endif
      /*
       * free memory
       */
      // free_discretization(&disc);
      // free_meshdata(&mesh);
    }
  }

  return out;
}
