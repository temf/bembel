// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "Visualize.hpp"

namespace Bembel {
namespace Vis {

double evaluate_laplace_rho_at(
    Bembel::Discretization<Bembel::LaplaceSingle> &myDisc,
    const Eigen::VectorXd &rho, int p, double x, double y) {
  Bembel::discretization &disc = myDisc.get_disc();
  const int ao = disc.a_o;
  const int elementNumber = disc.mesh->nf;
  const int bigSize = elementNumber * ao * ao;
  const int patchnum = disc.mesh->geom->size();
  const int M = disc.mesh->M;
  const int n = (1 << M);

  const double h = 1. / n;

  double *smallVec = Bembel::eigen2ptr(rho);
  double *longVec;
  longVec = (double *)calloc(bigSize, sizeof(double));

  proj_distr_et(&disc, smallVec, longVec);

  const Bembel::et_node *E = disc.mesh->E.patch[0];
  while (E->son[0]) E = E->son[0];

  double out = 0;
  for (int element = 0; element < n * n * patchnum; element++) {
    if (E[element].patch == p) {
      const double bot_x = (double)E[element].index_s / n;
      const double top_x = bot_x + h;
      if (bot_x <= x && top_x >= x) {
        const double bot_y = (double)E[element].index_t / n;
        const double top_y = bot_y + h;
        if (bot_y <= y && top_y >= y) {
          const double localized_x = Spl::rescale(x, bot_x, top_x);
          const double localized_y = Spl::rescale(y, bot_y, top_y);
          const int shift = element * ao * ao;
          double basis_vals[ao * ao];
          disc.phiphi(basis_vals,
                      Bembel::vector2_make(localized_x, localized_y));
          for (int i = 0; i < ao * ao; i++) {
            out += basis_vals[i] * longVec[shift + i];
          }
          free(longVec);
          free(smallVec);
          return out;
        }
      }
    }
  }
  std::cout << std::endl << "x " << x << "\ny " << y << std::endl;
  assert(false &&
         "If this happens, please write a bug report. This means that an "
         "evaluation point could not be assigned an element. However, this "
         "must always be the case.");

  return out;
}

std::complex<double> evaluate_helmholtz_rho_at(
    Bembel::Discretization<Bembel::HelmholtzSingle> &myDisc,
    const Eigen::VectorXcd &rho, int p, double x, double y) {
  Bembel::discretization &disc = myDisc.get_disc();
  const int ao = disc.a_o;
  const int elementNumber = disc.mesh->nf;
  const int bigSize = elementNumber * ao * ao;
  const int patchnum = disc.mesh->geom->size();
  const int M = disc.mesh->M;
  const int n = (1 << M);

  const double h = 1. / n;

  double *smallVec = Bembel::eigen2cmplxptr(rho);
  double *longVec;
  longVec = (double *)calloc(bigSize * 2, sizeof(double));

  proj_distr_et(&disc, smallVec, longVec);

  const Bembel::et_node *E = disc.mesh->E.patch[0];
  while (E->son[0]) E = E->son[0];

  double out_real = 0;
  double out_imag = 0;

  for (int element = 0; element < n * n * patchnum; element++) {
    if (E[element].patch == p) {
      const double bot_x = (double)E[element].index_s / n;
      const double top_x = bot_x + h;
      if (bot_x <= x && top_x > x) {
        const double bot_y = (double)E[element].index_t / n;
        const double top_y = bot_y + h;
        if (bot_y <= y && top_y > y) {
          const double localized_x = Spl::rescale(x, bot_x, top_x);
          const double localized_y = Spl::rescale(y, bot_y, top_y);
          const int shift = element * ao * ao;
          double basis_vals[ao * ao];
          disc.phiphi(basis_vals,
                      Bembel::vector2_make(localized_x, localized_y));
          for (int i = 0; i < ao * ao; i++) {
            out_real += basis_vals[i] * longVec[shift + i];
            out_imag += basis_vals[i] * longVec[bigSize + shift + i];
          }
          free(longVec);
          free(smallVec);
          return std::complex<double>(out_real, out_imag);
        }
      }
    }
  }
  assert(false &&
         "If this happens, please write a bug report. This means that an "
         "evaluation point could not be assigned an element. However, this "
         "must always be the case.");

  return std::complex<double>(0, 0);
}

Eigen::Vector3cd evaluate_maxwell_rho_at(
    Bembel::Discretization<Bembel::MaxwellSingle> &myDisc,
    const Eigen::VectorXcd &rho, int p, double x, double y) {
  Bembel::discretization &disc = myDisc.get_disc();
  std::vector<Spl::Patch> &geom = myDisc.get_plain_patchdata();
  const int ao = disc.a_o;
  const int elementNumber = disc.mesh->nf;
  const int bigSize = elementNumber * ao * ao;
  const int patchnum = disc.mesh->geom->size();
  const int M = disc.mesh->M;
  const int n = (1 << M);

  const double h = 1. / n;

  double *smallVec = Bembel::eigen2cmplxptr(rho);
  double *longVec;
  longVec = (double *)calloc(bigSize * 4, sizeof(double));

  proj_distr_et(&disc, smallVec, longVec);

  const Bembel::et_node *E = disc.mesh->E.patch[0];
  while (E->son[0]) E = E->son[0];

  double out_real_1 = 0;
  double out_imag_1 = 0;

  double out_real_2 = 0;
  double out_imag_2 = 0;

  for (int element = 0; element < n * n * patchnum; element++) {
    if (E[element].patch == p) {
      const double bot_x = (double)E[element].index_s / n;
      const double top_x = bot_x + h;
      if (bot_x <= x && top_x > x) {
        const double bot_y = (double)E[element].index_t / n;
        const double top_y = bot_y + h;
        if (bot_y <= y && top_y > y) {
          const double localized_x = Spl::rescale(x, bot_x, top_x);
          const double localized_y = Spl::rescale(y, bot_y, top_y);
          const int shift = element * ao * ao;
          double basis_vals[ao * ao];
          disc.phiphi(basis_vals,
                      Bembel::vector2_make(localized_x, localized_y));
          for (int i = 0; i < ao * ao; i++) {
            out_real_1 += basis_vals[i] * longVec[shift + i];
            out_imag_1 += basis_vals[i] * longVec[bigSize + shift + i];
            out_real_2 += basis_vals[i] * longVec[bigSize * 2 + shift + i];
            out_imag_2 += basis_vals[i] * longVec[bigSize * 3 + shift + i];
          }
          free(longVec);
          free(smallVec);
          Eigen::Matrix<double, 3, 2> jac = geom[p].jacobian(x, y);
          double surfmeas = geom[p].evaln(x, y).norm();
          return (jac.col(0) * std::complex<double>(out_real_1, out_imag_1) +
                 jac.col(1) * std::complex<double>(out_real_2, out_imag_2))/surfmeas;
        }
      }
    }
  }
  assert(false &&
         "If this happens, please write a bug report. This means that an "
         "evaluation point could not be assigned an element. However, this "
         "must always be the case.");
  Eigen::Matrix<double, 3, 2> jac = geom[p].jacobian(x, y);

  return jac.col(0) * std::complex<double>(out_real_1, out_imag_1) +
         jac.col(1) * std::complex<double>(out_real_2, out_imag_2);
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
/////    From here on follow the vtk-routines    /////
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

void plot_geometry_only(const std::vector<Spl::Patch> &geoms, int num,
                        const char *name) {
  using namespace Spl;
  check_geometry(geoms);
  num = (1 << num);
  num = num + 1;
  const double h = 1. / (num - 1);
  const int patchnum = geoms.size();
  const std::vector<double> grid = [&](int k) {
    std::vector<double> out;
    for (int i = 0; i < k; i++) {
      out.push_back(i * h);
    }
    out.push_back(1);
    out.shrink_to_fit();
    return out;
  }(num);
  const std::vector<double> center = [&](int k) {
    std::vector<double> out;
    for (int i = 0; i < k; i++) {
      out.push_back(i * h + .5 * h);
    }
    out.shrink_to_fit();
    return out;
  }(num);

  const int onedgridsize = grid.size();
  const int twodgridsize = onedgridsize * onedgridsize;
  const int onedcellsize = center.size();
  const int twodcellsize = onedcellsize * onedcellsize;

  std::vector<Eigen::Vector3d> pts;
  pts.reserve(twodgridsize * patchnum);

  std::vector<std::array<int, 4>> cells;
  std::vector<Eigen::Vector3d> normals;
  std::vector<int> patchId;
  std::vector<Eigen::Vector3d> dx, dy;

  cells.reserve(twodcellsize * patchnum);
  normals.reserve(twodcellsize * patchnum);
  patchId.reserve(twodcellsize * patchnum);
  dx.reserve(twodcellsize * patchnum);
  dy.reserve(twodcellsize * patchnum);

  for (int i = 0; i < patchnum; i++) {
    const Patch &g = geoms[i];
    for (int ix = 0; ix < onedgridsize; ix++)
      for (int iy = 0; iy < onedgridsize; iy++) {
        const auto x = grid[ix];
        const auto y = grid[iy];
        pts.push_back(g.eval(x, y));
      }
    for (int ix = 0; ix < onedcellsize; ix++)
      for (int iy = 0; iy < onedcellsize; iy++) {
        const auto x = center[ix];
        const auto y = center[iy];
        cells.push_back({
            iy + ix * onedgridsize + i * twodgridsize,            //
            iy + ix * onedgridsize + 1 + i * twodgridsize,        //
            iy + (ix + 1) * onedgridsize + 1 + i * twodgridsize,  //
            iy + (ix + 1) * onedgridsize + i * twodgridsize       //
        });                                                       //
        patchId.push_back(i);
        auto jac = g.jacobian(x, y);
        normals.push_back(g.evaln(x, y));  // jac-core does not see the flip
        dx.push_back(jac.col(0).eval());
        dy.push_back(jac.col(1).eval());
      }
  }

  const int cz = cells.size();
  const int gz = pts.size();

  // std::cout << "                           patches: " << geoms.size() <<
  // "\n"; std::cout << "cz   " << cz << "\n"; std::cout << "gz   " << gz <<
  // "\n";

  FILE *g = NULL;
  g = fopen(name, "w");

  /*
   * vtk format header
   */
  fprintf(g, "# vtk DataFile Version 3.1\n");
  fprintf(g, "this file hopefully represents my surface now\n");
  fprintf(g, "ASCII\n");
  fprintf(g, "DATASET UNSTRUCTURED_GRID\n");

  /*
   * print point list
   */
  fprintf(g, "POINTS %d FLOAT\n", gz);
  for (int i = 0; i < gz; ++i)
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", pts[i](0), pts[i](1), pts[i](2));
  fprintf(g, "\n");

  /*
   * print element list
   */
  fprintf(g, "CELLS %d %d\n", cz, 5 * cz);
  for (int i = 0; i < cz; ++i)
    fprintf(g, "%d %d %d %d %d\n", 4, cells[i][0], cells[i][1], cells[i][2],
            cells[i][3]);
  fprintf(g, "\n");

  fprintf(g, "CELL_TYPES %d\n", cz);
  for (int i = 0; i < cz; ++i) fprintf(g, "%d\n", 9);
  fprintf(g, "\n");

  /*
   * print z-values of the geometry and solved density for visualization
   */
  fprintf(g, "CELL_DATA %d\n", cz);
  fprintf(g, "VECTORS normals FLOAT\n");
  for (int i = 0; i < cz; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", normals[i](0), normals[i](1),
            normals[i](2));
  }

  fprintf(g, "\n");
  fprintf(g, "VECTORS ds FLOAT\n");
  for (int i = 0; i < cz; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", dx[i](0), dx[i](1), dx[i](2));
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS dt FLOAT\n");
  for (int i = 0; i < cz; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", dy[i](0), dy[i](1), dy[i](2));
  }

  fprintf(g, "SCALARS patch FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (int i = 0; i < cz; ++i) fprintf(g, "%d\n", patchId[i]);
  fprintf(g, "\n");
  fclose(g);

  return;
}

void plot_rho_laplace(Bembel::Discretization<Bembel::LaplaceSingle> &myDisc,
                      const Eigen::VectorXd &rho, int num, const char *name) {
  using namespace Spl;
  std::vector<Spl::Patch> &geoms = myDisc.get_plain_patchdata();

  check_geometry(geoms);
  num = (1 << num);
  num = num + 1;
  const double h = 1. / (num - 1);
  const int patchnum = geoms.size();
  const std::vector<double> grid = [&](int k) {
    std::vector<double> out;
    for (int i = 0; i < k; i++) {
      out.push_back(i * h);
      std::cout << out[i] << std::endl;
    }
    out.push_back(1);
    out.shrink_to_fit();
    return out;
  }(num);
  const std::vector<double> center = [&](int k) {
    std::vector<double> out;
    for (int i = 0; i < k - 1; i++) {
      out.push_back(i * h + .5 * h);
    }
    out.shrink_to_fit();
    return out;
  }(num);

  const int onedgridsize = grid.size();
  const int twodgridsize = onedgridsize * onedgridsize;
  const int onedcellsize = center.size();
  const int twodcellsize = onedcellsize * onedcellsize;

  std::vector<Eigen::Vector3d> pts;
  pts.reserve(twodgridsize * patchnum);

  std::vector<std::array<int, 4>> cells;
  std::vector<Eigen::Vector3d> normals;
  std::vector<int> patchId;
  std::vector<Eigen::Vector3d> dx, dy;
  std::vector<double> rho_pointwise;

  cells.reserve(twodcellsize * patchnum);
  normals.reserve(twodcellsize * patchnum);
  patchId.reserve(twodcellsize * patchnum);
  dx.reserve(twodcellsize * patchnum);
  dy.reserve(twodcellsize * patchnum);
  rho_pointwise.reserve(twodcellsize * patchnum);

  for (int i = 0; i < patchnum; i++) {
    const Patch &g = geoms[i];
    for (int ix = 0; ix < onedgridsize; ix++)
      for (int iy = 0; iy < onedgridsize; iy++) {
        const auto x = grid[ix];
        const auto y = grid[iy];
        pts.push_back(g.eval(x, y));
      }
    for (int ix = 0; ix < onedcellsize; ix++)
      for (int iy = 0; iy < onedcellsize; iy++) {
        const auto x = center[ix];
        const auto y = center[iy];
        cells.push_back({
            iy + ix * onedgridsize + i * twodgridsize,            //
            iy + ix * onedgridsize + 1 + i * twodgridsize,        //
            iy + (ix + 1) * onedgridsize + 1 + i * twodgridsize,  //
            iy + (ix + 1) * onedgridsize + i * twodgridsize       //
        });                                                       //
        patchId.push_back(i);
        auto jac = g.jacobian(x, y);
        normals.push_back(g.evaln(x, y));  // jac-core does not see the flip
        dx.push_back(jac.col(0).eval());
        dy.push_back(jac.col(1).eval());

        rho_pointwise.push_back(evaluate_laplace_rho_at(myDisc, rho, i, x, y));
      }
  }

  const int cz = cells.size();
  const int gz = pts.size();

  // std::cout << "                           patches: " << geoms.size() <<
  // "\n"; std::cout << "cz   " << cz << "\n"; std::cout << "gz   " << gz <<
  // "\n";

  FILE *g = NULL;
  g = fopen(name, "w");

  /*
   * vtk format header
   */
  fprintf(g, "# vtk DataFile Version 3.1\n");
  fprintf(g, "this file hopefully represents my surface now\n");
  fprintf(g, "ASCII\n");
  fprintf(g, "DATASET UNSTRUCTURED_GRID\n");

  /*
   * print point list
   */
  fprintf(g, "POINTS %d FLOAT\n", gz);
  for (int i = 0; i < gz; ++i)
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", pts[i](0), pts[i](1), pts[i](2));
  fprintf(g, "\n");

  /*
   * print element list
   */
  fprintf(g, "CELLS %d %d\n", cz, 5 * cz);
  for (int i = 0; i < cz; ++i)
    fprintf(g, "%d %d %d %d %d\n", 4, cells[i][0], cells[i][1], cells[i][2],
            cells[i][3]);
  fprintf(g, "\n");

  fprintf(g, "CELL_TYPES %d\n", cz);
  for (int i = 0; i < cz; ++i) fprintf(g, "%d\n", 9);
  fprintf(g, "\n");

  /*
   * print z-values of the geometry and solved density for visualization
   */
  fprintf(g, "CELL_DATA %d\n", cz);
  fprintf(g, "VECTORS normals FLOAT\n");
  for (int i = 0; i < cz; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", normals[i](0), normals[i](1),
            normals[i](2));
  }

  fprintf(g, "\n");
  fprintf(g, "VECTORS ds FLOAT\n");
  for (int i = 0; i < cz; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", dx[i](0), dx[i](1), dx[i](2));
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS dt FLOAT\n");
  for (int i = 0; i < cz; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", dy[i](0), dy[i](1), dy[i](2));
  }

  fprintf(g, "SCALARS patch FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (int i = 0; i < cz; ++i) fprintf(g, "%d\n", patchId[i]);
  fprintf(g, "\n");

  fprintf(g, "SCALARS rho FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (int i = 0; i < cz; ++i) fprintf(g, "%20.16f\n", rho_pointwise[i]);
  fprintf(g, "\n");

  fclose(g);

  return;
}

///////////////////////////////////////////////////////////////////////////////

void plot_rho_helmholtz(Bembel::Discretization<Bembel::HelmholtzSingle> &myDisc,
                        const Eigen::VectorXcd &rho, int num,
                        const char *name) {
  using namespace Spl;
  std::vector<Spl::Patch> &geoms = myDisc.get_plain_patchdata();

  check_geometry(geoms);
  num = (1 << num);
  num = num + 1;
  const double h = 1. / (num - 1);
  const int patchnum = geoms.size();
  const std::vector<double> grid = [&](int k) {
    std::vector<double> out;
    for (int i = 0; i < num; i++) {
      out.push_back(i * h);
    }
    out.push_back(1);
    out.shrink_to_fit();
    return out;
  }(num);
  const std::vector<double> center = [&](int k) {
    std::vector<double> out;
    for (int i = 0; i < k - 1; i++) {
      out.push_back(i * h + .5 * h);
    }
    out.shrink_to_fit();
    return out;
  }(num);

  const int onedgridsize = grid.size();
  const int twodgridsize = onedgridsize * onedgridsize;
  const int onedcellsize = center.size();
  const int twodcellsize = onedcellsize * onedcellsize;

  std::vector<Eigen::Vector3d> pts;
  pts.reserve(twodgridsize * patchnum);

  std::vector<std::array<int, 4>> cells;
  std::vector<Eigen::Vector3d> normals;
  std::vector<int> patchId;
  std::vector<Eigen::Vector3d> dx, dy;
  std::vector<std::complex<double>> rho_pointwise;

  cells.reserve(twodcellsize * patchnum);
  normals.reserve(twodcellsize * patchnum);
  patchId.reserve(twodcellsize * patchnum);
  dx.reserve(twodcellsize * patchnum);
  dy.reserve(twodcellsize * patchnum);
  rho_pointwise.reserve(twodcellsize * patchnum);

  for (int i = 0; i < patchnum; i++) {
    const Patch &g = geoms[i];
    for (int ix = 0; ix < onedgridsize; ix++)
      for (int iy = 0; iy < onedgridsize; iy++) {
        const auto x = grid[ix];
        const auto y = grid[iy];
        pts.push_back(g.eval(x, y));
      }
    for (int ix = 0; ix < onedcellsize; ix++)
      for (int iy = 0; iy < onedcellsize; iy++) {
        const auto x = center[ix];
        const auto y = center[iy];
        cells.push_back({
            iy + ix * onedgridsize + i * twodgridsize,            //
            iy + ix * onedgridsize + 1 + i * twodgridsize,        //
            iy + (ix + 1) * onedgridsize + 1 + i * twodgridsize,  //
            iy + (ix + 1) * onedgridsize + i * twodgridsize       //
        });                                                       //
        patchId.push_back(i);
        auto jac = g.jacobian(x, y);
        normals.push_back(g.evaln(x, y));  // jac-core does not see the flip
        dx.push_back(jac.col(0).eval());
        dy.push_back(jac.col(1).eval());
        rho_pointwise.push_back(
            evaluate_helmholtz_rho_at(myDisc, rho, i, x, y));
      }
  }

  const int cz = cells.size();
  const int gz = pts.size();

  // std::cout << "                           patches: " << geoms.size() <<
  // "\n"; std::cout << "cz   " << cz << "\n"; std::cout << "gz   " << gz <<
  // "\n";

  FILE *g = NULL;
  g = fopen(name, "w");

  /*
   * vtk format header
   */
  fprintf(g, "# vtk DataFile Version 3.1\n");
  fprintf(g, "this file hopefully represents my surface now\n");
  fprintf(g, "ASCII\n");
  fprintf(g, "DATASET UNSTRUCTURED_GRID\n");

  /*
   * print point list
   */
  fprintf(g, "POINTS %d FLOAT\n", gz);
  for (int i = 0; i < gz; ++i)
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", pts[i](0), pts[i](1), pts[i](2));
  fprintf(g, "\n");

  /*
   * print element list
   */
  fprintf(g, "CELLS %d %d\n", cz, 5 * cz);
  for (int i = 0; i < cz; ++i)
    fprintf(g, "%d %d %d %d %d\n", 4, cells[i][0], cells[i][1], cells[i][2],
            cells[i][3]);
  fprintf(g, "\n");

  fprintf(g, "CELL_TYPES %d\n", cz);
  for (int i = 0; i < cz; ++i) fprintf(g, "%d\n", 9);
  fprintf(g, "\n");

  /*
   * print z-values of the geometry and solved density for visualization
   */
  fprintf(g, "CELL_DATA %d\n", cz);
  fprintf(g, "VECTORS normals FLOAT\n");
  for (int i = 0; i < cz; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", normals[i](0), normals[i](1),
            normals[i](2));
  }

  fprintf(g, "\n");
  fprintf(g, "VECTORS ds FLOAT\n");
  for (int i = 0; i < cz; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", dx[i](0), dx[i](1), dx[i](2));
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS dt FLOAT\n");
  for (int i = 0; i < cz; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", dy[i](0), dy[i](1), dy[i](2));
  }

  fprintf(g, "SCALARS patch FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (int i = 0; i < cz; ++i) fprintf(g, "%d\n", patchId[i]);
  fprintf(g, "\n");

  fprintf(g, "SCALARS rho_real FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (int i = 0; i < cz; ++i) fprintf(g, "%20.16f\n", rho_pointwise[i].real());
  fprintf(g, "\n");

  fprintf(g, "SCALARS rho_imag FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (int i = 0; i < cz; ++i) fprintf(g, "%20.16f\n", rho_pointwise[i].imag());
  fprintf(g, "\n");

  fclose(g);

  return;
}

//////////////////////////////////////////////////////////////////

void plot_rho_maxwell(Bembel::Discretization<Bembel::MaxwellSingle> &myDisc,
                      const Eigen::VectorXcd &rho, int num, const char *name) {
  using namespace Spl;
  std::vector<Spl::Patch> &geoms = myDisc.get_plain_patchdata();

  check_geometry(geoms);
  num = (1 << num);
  num = num + 1;
  const double h = 1. / (num - 1);
  const int patchnum = geoms.size();
  const std::vector<double> grid = [&](int k) {
    std::vector<double> out;
    for (int i = 0; i < k; i++) {
      out.push_back(i * h);
    }
    out.push_back(1);
    out.shrink_to_fit();
    return out;
  }(num);
  const std::vector<double> center = [&](int k) {
    std::vector<double> out;
    for (int i = 0; i < k - 1; i++) {
      out.push_back(i * h + .5 * h);
    }
    out.shrink_to_fit();
    return out;
  }(num);

  const int onedgridsize = grid.size();
  const int twodgridsize = onedgridsize * onedgridsize;
  const int onedcellsize = center.size();
  const int twodcellsize = onedcellsize * onedcellsize;

  std::vector<Eigen::Vector3d> pts;
  pts.reserve(twodgridsize * patchnum);

  std::vector<std::array<int, 4>> cells;
  std::vector<Eigen::Vector3d> normals;
  std::vector<int> patchId;
  std::vector<Eigen::Vector3d> dx, dy;
  std::vector<Eigen::Vector3cd> rho_pointwise;

  cells.reserve(twodcellsize * patchnum);
  normals.reserve(twodcellsize * patchnum);
  patchId.reserve(twodcellsize * patchnum);
  dx.reserve(twodcellsize * patchnum);
  dy.reserve(twodcellsize * patchnum);
  rho_pointwise.reserve(twodcellsize * patchnum);

  for (int i = 0; i < patchnum; i++) {
    const Patch &g = geoms[i];
    for (int ix = 0; ix < onedgridsize; ix++)
      for (int iy = 0; iy < onedgridsize; iy++) {
        const auto x = grid[ix];
        const auto y = grid[iy];
        pts.push_back(g.eval(x, y));
      }
    for (int ix = 0; ix < onedcellsize; ix++)
      for (int iy = 0; iy < onedcellsize; iy++) {
        const auto x = center[ix];
        const auto y = center[iy];
        cells.push_back({
            iy + ix * onedgridsize + i * twodgridsize,            //
            iy + ix * onedgridsize + 1 + i * twodgridsize,        //
            iy + (ix + 1) * onedgridsize + 1 + i * twodgridsize,  //
            iy + (ix + 1) * onedgridsize + i * twodgridsize       //
        });                                                       //
        patchId.push_back(i);
        auto jac = g.jacobian(x, y);
        normals.push_back(g.evaln(x, y));  // jac-core does not see the
        dx.push_back(jac.col(0).eval());
        dy.push_back(jac.col(1).eval());
        rho_pointwise.push_back(evaluate_maxwell_rho_at(myDisc, rho, i, x, y));
      }
  }

  const int cz = cells.size();
  const int gz = pts.size();

  // std::cout << "                           patches: " << geoms.size() <<
  // "\n"; std::cout << "cz   " << cz << "\n"; std::cout << "gz   " << gz <<
  // "\n";

  FILE *g = NULL;
  g = fopen(name, "w");

  /*
   * vtk format header
   */
  fprintf(g, "# vtk DataFile Version 3.1\n");
  fprintf(g, "this file hopefully represents my surface now\n");
  fprintf(g, "ASCII\n");
  fprintf(g, "DATASET UNSTRUCTURED_GRID\n");

  /*
   * print point list
   */
  fprintf(g, "POINTS %d FLOAT\n", gz);
  for (int i = 0; i < gz; ++i)
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", pts[i](0), pts[i](1), pts[i](2));
  fprintf(g, "\n");

  /*
   * print element list
   */
  fprintf(g, "CELLS %d %d\n", cz, 5 * cz);
  for (int i = 0; i < cz; ++i)
    fprintf(g, "%d %d %d %d %d\n", 4, cells[i][0], cells[i][1], cells[i][2],
            cells[i][3]);
  fprintf(g, "\n");

  fprintf(g, "CELL_TYPES %d\n", cz);
  for (int i = 0; i < cz; ++i) fprintf(g, "%d\n", 9);
  fprintf(g, "\n");

  /*
   * print z-values of the geometry and solved density for visualization
   */
  fprintf(g, "CELL_DATA %d\n", cz);
  fprintf(g, "VECTORS normals FLOAT\n");
  for (int i = 0; i < cz; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", normals[i](0), normals[i](1),
            normals[i](2));
  }

  fprintf(g, "\n");
  fprintf(g, "VECTORS ds FLOAT\n");
  for (int i = 0; i < cz; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", dx[i](0), dx[i](1), dx[i](2));
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS dt FLOAT\n");
  for (int i = 0; i < cz; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", dy[i](0), dy[i](1), dy[i](2));
  }

  fprintf(g, "SCALARS patch FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (int i = 0; i < cz; ++i) fprintf(g, "%d\n", patchId[i]);
  fprintf(g, "\n");

  fprintf(g, "VECTORS rho_real FLOAT\n");
  for (int i = 0; i < cz; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", rho_pointwise[i](0).real(),
            rho_pointwise[i](1).real(), rho_pointwise[i](2).real());
  }
  fprintf(g, "\n");

  fprintf(g, "VECTORS rho_imag FLOAT\n");
  for (int i = 0; i < cz; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", rho_pointwise[i](0).imag(),
            rho_pointwise[i](1).imag(), rho_pointwise[i](2).imag());
  }
  fprintf(g, "\n");

  fclose(g);

  return;
}
}  // namespace Vis
}  // namespace Bembel
