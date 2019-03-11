// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "discretization.h"
#include "glue.h"
#include "spline/basis.h"
#include "spline/deBoor.h"
#include "spline/knots.h"
#include "spline/localise.h"

// #ifdef _DBG_PROJECTOR_
#include <fstream>
// #endif

namespace Bembel {

int init_projector_Laplace(discretization *disc, et_node *E, const int M,
                           const int pp1, const int kntrep) {
  const auto info = init_projector_base(disc, pp1, pp1, kntrep);

  const int size = info.vals.size();
  std::vector<Eigen::Triplet<double>> trips;

  for (int k = 0; k < size; k++) {
    trips.push_back(
        Eigen::Triplet<double>(info.rows[k], info.cols[k], info.vals[k]));
  }

  Eigen::SparseMatrix<double> Proj(info.num_dc_dof, info.num_c_dof);
  Proj.setFromTriplets(trips.begin(), trips.end());

  sparse **pT = &disc->T;
  sparse *T;

  T = *pT = (sparse *)calloc(1, sizeof(sparse));

  init_sparse(T, info.num_dc_dof, info.num_c_dof, pp1 * pp1);

  for (int j = 0; j < size; j++)
    set_sparse<double>(T, info.rows[j], info.cols[j], info.vals[j]);

  // init_sparse(T, info.num_dc_dof, info.num_dc_dof, 1);
  // for (int j = 0; j < info.num_dc_dof; j++)
  // set_sparse<double>(T, j, j, 1);
  // return info.num_dc_dof;

  return info.num_c_dof;
}

int init_projector_Laplace_cont(discretization *disc, et_node *E, const int M,
                                const int pp1, const int kntrep) {
  const auto info = init_projector_base(disc, pp1, pp1, kntrep);

  const int size = info.vals.size();
  std::vector<Eigen::Triplet<double>> trips;

  for (int k = 0; k < size; k++) {
    trips.push_back(
        Eigen::Triplet<double>(info.rows[k], info.cols[k], info.vals[k]));
  }

  Eigen::SparseMatrix<double> noGlue(info.num_dc_dof, info.num_c_dof);
  noGlue.setFromTriplets(trips.begin(), trips.end());
  Eigen::SparseMatrix<double> Glue = make_glue_matrix_laplace(
      disc, pp1, kntrep, M,
      info.num_c_dof);  // I think the pp1 input is no longer used.
  Eigen::SparseMatrix<double> project = (noGlue * Glue.transpose()).eval();
  // Eigen::SparseMatrix<double> Proj = (noGlue).eval();

  sparse **pT = &disc->T;
  sparse *T;

  T = *pT = (sparse *)calloc(1, sizeof(sparse));
  init_sparse(T, project.rows(), project.cols(), 2 * (pp1) * (pp1));

  for (int k = 0; k < project.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(project, k); it; ++it) {
      set_sparse<double>(T, it.row(), it.col(), it.value());
    }

  return project.cols();
}

int init_projector_Helmholtz(discretization *disc, et_node *E, const int M,
                             const int pp1, const int kntrep) {
  const auto info = init_projector_base(disc, pp1, pp1, kntrep);

  const int size = info.vals.size();

  const int num_c_dof = info.num_c_dof;
  const int num_dc_dof = info.num_dc_dof;

  sparse **pT = &disc->T;
  sparse *T;

  T = *pT = (sparse *)calloc(1, sizeof(sparse));
  init_sparse(T, info.num_dc_dof * 2, info.num_c_dof * 2, pp1 * pp1);

  for (int j = 0; j < size; j++)
    set_sparse<double>(T, info.rows[j], info.cols[j], info.vals[j]);
  for (int j = 0; j < size; j++)
    set_sparse<double>(T, info.rows[j] + num_dc_dof, info.cols[j] + num_c_dof,
                       info.vals[j]);

  return info.num_c_dof + info.num_c_dof;  // ????
}

int init_projector_Maxwell(discretization *disc, et_node *E, const int M,
                           const int pp1, int kntrep) {
  // std::cout << "KNTREP IS "<< kntrep << "\n";

  // kntrep should be between 1 and pp1-1
  kntrep = std::max(std::min(kntrep, pp1 - 1), 1);

  auto info = init_projector_base(disc, pp1, pp1 - 1, kntrep);

  const int size = info.vals.size();

  const int num_c_dof = info.num_c_dof;
  const int num_dc_dof = info.num_dc_dof;

  sparse **pT = &disc->T;
  sparse *T;

  std::vector<Eigen::Triplet<double>> trips;

  for (int k = 0; k < size; k++) {
    trips.push_back(Eigen::Triplet<double>(info.cols[k] + 0 * num_c_dof,
                                           info.rows[k] + 0 * num_dc_dof,
                                           info.vals[k]));
    trips.push_back(Eigen::Triplet<double>(info.cols[k] + 2 * num_c_dof,
                                           info.rows[k] + 1 * num_dc_dof,
                                           info.vals[k]));
  }

  info = init_projector_base(disc, pp1 - 1, pp1, kntrep);

  assert(info.num_c_dof == num_c_dof and info.num_dc_dof == num_dc_dof);

  for (int k = 0; k < size; k++) {
    trips.push_back(Eigen::Triplet<double>(info.cols[k] + 1 * num_c_dof,
                                           info.rows[k] + 2 * num_dc_dof,
                                           info.vals[k]));
    trips.push_back(Eigen::Triplet<double>(info.cols[k] + 3 * num_c_dof,
                                           info.rows[k] + 3 * num_dc_dof,
                                           info.vals[k]));
  }

  Eigen::SparseMatrix<double> NoGlue(num_c_dof * 4, num_dc_dof * 4);
  NoGlue.setFromTriplets(trips.begin(), trips.end());

  Eigen::SparseMatrix<double> Glue = make_glue_matrix_maxwell(
      disc, pp1, kntrep, M,
      num_c_dof * 4);  // I think the pp1 input is no longer used.

#ifdef _BEMBEL_PRINT_INFO_
  std::cout << "                     Dofs: Final Update\n";
  std::cout << "                           Big Matrix = " << NoGlue.cols()
            << "\n                     Small Matrix (not glued) = "
            << NoGlue.rows()
            << "\n                     Small Matrix (glued) = " << Glue.rows()
            << "\n";
  std::cout << "                           knotrepetition = " << kntrep << "\n";
#endif
  Eigen::SparseMatrix<double> project = (Glue * NoGlue).transpose();
  // Eigen::SparseMatrix<double> project = (NoGlue).transpose();

  T = *pT = (sparse *)calloc(1, sizeof(sparse));
  init_sparse(T, project.rows(), project.cols(), 2 * (pp1) * (pp1));

  for (int k = 0; k < project.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(project, k); it; ++it) {
      set_sparse<double>(T, it.row(), it.col(), it.value());
    }

    // std::cout << "                           " << project.cols()
    // << " dofs after projection assembly\n";

#if 0
  std::ofstream dbg1;
  dbg1.open("glue.dbg.log");
  dbg1 << Eigen::MatrixXd(Glue);
  dbg1.close();

  std::ofstream dbg2;
  dbg2.open("noglue.dbg.log");
  dbg2 << Eigen::MatrixXd(NoGlue);
  dbg2.close();

  std::ofstream dbg3;
  dbg3.open("project.dbg.log");
  dbg3 << Eigen::MatrixXd(project);
  dbg3.close();
#endif

  return project.cols();
}

/**
 *  \brief         Initializes the projector matrix used in proj_distr_et()
 *                 and proj_restr_et(). Does so by interpolation of the TP
 *           B-Spline space.
 *
 *  \param[in]     E              Element list
 *  \param[in]     M              Level
 *
 *  \return     na             Number of ansatz functions
 *  \attention Here we assume that our de Rham sequence starts with a space of
 * the same degree in every TP-direction. \author        Felix Wolf
 */
_proj_info init_projector_base(discretization *disc, const int pp1x,
                               const int pp1y, const int kntrep_in) {
  using namespace Spl;

  const et_node *E = disc->mesh->E.patch[0];
  const int M = disc->mesh->M;
  const int maxp = std::max(pp1x, pp1y);
  const int minp = std::min(pp1x, pp1y);
  const int kntrep = (maxp <= kntrep_in) ? minp : kntrep_in;
  const int n = (1 << M);
  const int patchnum = disc->mesh->geom->size();

  // At the moment, we allow a difference of one only, i.e., we start our de
  // Rham sequence with a space of the same degree in every TP-direction.
  assert(std::abs(maxp - minp) <= 1);

  // Now, we construct the continuous spaces which will be the preimage of the
  // projector. n-1 is passed, since n = number_elements but n-1 = number of
  // interior knots.
  std::vector<double> c_space_knot_x = make_unif_knots(pp1x, n - 1, kntrep);
  std::vector<double> c_space_knot_y = make_unif_knots(pp1y, n - 1, kntrep);
  const int c_space_dim_x = c_space_knot_x.size() - pp1x;
  const int c_space_dim_y = c_space_knot_y.size() - pp1y;
  const int c_space_dim = c_space_dim_x * c_space_dim_y;

  _proj_info out;

  out.cols.reserve(c_space_dim * maxp * maxp * patchnum);
  out.rows.reserve(c_space_dim * maxp * maxp * patchnum);
  out.vals.reserve(c_space_dim * maxp * maxp * patchnum);

  std::vector<double> mask = make_interpolation_mask(maxp);
  const int masksize = mask.size();

  // Now we assemble our system which we use for the interpolation
  assert(masksize == maxp && "projector.cpp: System needs to be square");
  Eigen::Matrix<double, -1, -1> system(masksize * masksize, maxp * maxp);

  // Here, we suddenly use the degree for basisevaluation, i.e., maxp-1. This is
  // confusing, but correct and tested.
  for (int iy = 0; iy < maxp; iy++) {
    std::vector<double> vals_y = evalBrnstnBasis(maxp - 1, mask[iy]);
    for (int ix = 0; ix < maxp; ix++) {
      std::vector<double> vals_x = evalBrnstnBasis(maxp - 1, mask[ix]);
      for (int jy = 0; jy < maxp; jy++) {
        for (int jx = 0; jx < maxp; jx++) {
          system(iy * masksize + ix, jy * maxp + jx) = vals_x[jx] * vals_y[jy];
        }
      }
    }
  }

  // Do LU-Decomposition of system
  const Eigen::PartialPivLU<Eigen::Matrix<double, -1, -1>> pplu(system);

  // find lowest level in element tree
  while (E->son[0]) E = E->son[0];

  for (int element = 0; element < n * n * patchnum; element++) {
    // s = x
    const double pos_x = (double)E[element].index_s / n;
    const double pos_y = (double)E[element].index_t / n;

    const double h = 1. / n;
    const double mid_x = pos_x + (.5 * h);
    const double mid_y = pos_y + (.5 * h);

    // The following block of code computes the indeces of all continuous basis
    // functions which have a support on the element of index 'element'
    std::vector<double> nonzero_dofs_x;
    std::vector<double> nonzero_dofs_y;
    nonzero_dofs_x.reserve(pp1x);
    nonzero_dofs_y.reserve(pp1y);
    for (int dof_y = 0; dof_y < c_space_dim_y; dof_y++) {
      if ((c_space_knot_y[dof_y] <= mid_y) and
          (c_space_knot_y[dof_y + pp1y] >= mid_y)) {
        nonzero_dofs_y.push_back(dof_y);
      }
    }
    for (int dof_x = 0; dof_x < c_space_dim_x; dof_x++) {
      if ((c_space_knot_x[dof_x] <= mid_x) and
          (c_space_knot_x[dof_x + pp1x] >= mid_x)) {
        nonzero_dofs_x.push_back(dof_x);
      }
    }
    nonzero_dofs_x.shrink_to_fit();
    nonzero_dofs_y.shrink_to_fit();

    // We now loop over all basis functions wich are non-zero on the element
    for (auto dof_y : nonzero_dofs_y) {
      for (auto dof_x : nonzero_dofs_x) {
        // The next few lines are to check if the element is within the support
        // of the basis function given by dof_x and dof_y.
        std::vector<double> local_mask_x(masksize);
        std::vector<double> local_mask_y(masksize);
        for (int i = 0; i < masksize; i++) {
          local_mask_x[i] = pos_x + (1.0 / n) * mask[i];
          local_mask_y[i] = pos_y + (1.0 / n) * mask[i];
          assert(local_mask_x[i] < 1 and local_mask_x[i] > 0);
          assert(local_mask_y[i] < 1 and local_mask_y[i] > 0);
        }

        // build unit vectors
        std::vector<double> c_coefs_x(c_space_dim_x, 0.0);
        std::vector<double> c_coefs_y(c_space_dim_y, 0.0);
        c_coefs_x[dof_x] = 1.;
        c_coefs_y[dof_y] = 1.;

        // evaluate rhs for interpolation problem
        std::vector<double> vals_x =
            deBoor(c_coefs_x, c_space_knot_x, local_mask_x);
        std::vector<double> vals_y =
            deBoor(c_coefs_y, c_space_knot_y, local_mask_y);

        // assemble and solve interpolation problem
        Eigen::VectorXd rhs(masksize * masksize);
        for (int iy = 0; iy < masksize; iy++)
          for (int ix = 0; ix < masksize; ix++)
            rhs(iy * masksize + ix) = vals_x[ix] * vals_y[iy];
        Eigen::VectorXd sol = pplu.solve(rhs);

        // insert entries into triplet
        for (int i = 0; i < maxp * maxp; i++) {
          if (std::abs(sol(i)) > 1e-7) {
            out.cols.push_back(E[element].patch * (c_space_dim) +
                               (c_space_dim_x * dof_y) + dof_x);
            out.rows.push_back((element * maxp * maxp) + i);
            out.vals.push_back(sol(i));
          }
        }
      }
    }
  }

  out.num_c_dof = c_space_dim_x * c_space_dim_y * patchnum;
  out.num_dc_dof = maxp * maxp * n * n * patchnum;

  out.cols.shrink_to_fit();
  out.rows.shrink_to_fit();
  out.vals.shrink_to_fit();

  assert((out.vals.size() == out.rows.size()) &&
         (out.rows.size() == out.cols.size()) &&
         "projector.cpp: you made a mistake, try again.");

  return out;
}
}  // namespace Bembel