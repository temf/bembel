// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "H_print_geometry.h"

namespace Bembel {
int print_geometry_Maxwell(discretization *disc, double *rho, char *dname,
                           double kappa[2]) {
  et_node *E = disc->mesh->E.patch[0];
  vector3 *P = disc->mesh->P;
  int n = disc->mesh->n;
  const int nf = disc->mesh->nf;
  // std::cout << "nf = "<< nf << "\n";
  const int np = disc->mesh->np;
  const int a_bs = disc->a_bs;
  int i = 0;
  // int np = 0; /* number of edge points, starting point */
  double h = 1. / n;
  double *surfmeas;
  vector2 s; /* left lower corner on square for zi */
  vector3 vec;
  vector3 vecc;
  vector3 vec2[2];
  vector3 *fmid;
  vector3 *normals;
  vector3 *dx;
  vector3 *dy;
  // parametrix *Chi = disc->mesh->geom->Chi;
  const geometry &Chi = *disc->mesh->geom;
  /*
   * for edge points
   */
  FILE *g = NULL;

  /*
   * get pointer to the leafs of the element tree
   */
  while (E->son[0]) E = E->son[0];

  fmid = (vector3 *)calloc(nf, sizeof(vector3));
  dx = (vector3 *)calloc(nf, sizeof(vector3));
  dy = (vector3 *)calloc(nf, sizeof(vector3));
  normals = (vector3 *)calloc(nf, sizeof(vector3));
  surfmeas = (double *)calloc(nf, sizeof(double));

  for (i = 0; i < nf; ++i) {
    s.x = E[i].index_s * h + 0.5 * h;
    s.y = E[i].index_t * h + 0.5 * h;
    normals[i] = Chi[E[i].patch].n_f(s);
    dx[i] = Chi[E[i].patch].df_dx(s);
    dy[i] = Chi[E[i].patch].df_dy(s);
    fmid[i] = Chi[E[i].patch].f(s);
  }

  for (i = 0; i < nf; ++i) {
    surfmeas[i] = vector3_norm(normals[i]);
    normals[i] = vector3_Smul(1. / surfmeas[i], normals[i]);
    dx[i] = vector3_Smul(1. / surfmeas[i], dx[i]);
    dy[i] = vector3_Smul(1. / surfmeas[i], dy[i]);
  }

  // Expand Rho
  double *bigrho;
  bigrho = (double *)calloc(nf * a_bs * 4, sizeof(double));
  proj_distr_et(disc, rho, bigrho);
  double rhovals[4 * nf];
  double phiphivec[a_bs];
  memset(phiphivec, 0, a_bs * sizeof(double));
  disc->phiphi(phiphivec, vector2_make(0.5, 0.5));
  for (i = 0; i < 4 * nf; i++) {
    rhovals[i] = myddot(a_bs, phiphivec, &bigrho[i * a_bs]);
  }
  double rhodvals[4 * nf];
  memset(phiphivec, 0, a_bs * sizeof(double));
  disc->phiphi_dx(phiphivec, vector2_make(0.5, 0.5));
  for (i = 0; i < 2 * nf; i++) {
    rhodvals[i] = myddot(a_bs, phiphivec, &bigrho[i * a_bs]);
  }
  memset(phiphivec, 0, a_bs * sizeof(double));
  disc->phiphi_dy(phiphivec, vector2_make(0.5, 0.5));
  for (i = 2 * nf; i < 4 * nf; i++) {
    rhodvals[i] = myddot(a_bs, phiphivec, &bigrho[i * a_bs]);
  }

  /*
   * visualize the surface of the geometry in vtk format
   */
  g = fopen(dname, "w");

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
  fprintf(g, "POINTS %d FLOAT\n", np);
  for (i = 0; i < np; ++i)
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", P[i].x, P[i].y, P[i].z);
  fprintf(g, "\n");

  /*
   * print element list
   */
  fprintf(g, "CELLS %d %d\n", nf, 5 * nf);
  for (i = 0; i < nf; ++i)
    fprintf(g, "%d %d %d %d %d\n", 4, E[i].vertex[0], E[i].vertex[1],
            E[i].vertex[2], E[i].vertex[3]);
  fprintf(g, "\n");

  fprintf(g, "CELL_TYPES %d\n", nf);
  for (i = 0; i < nf; ++i) fprintf(g, "%d\n", 9);
  fprintf(g, "\n");

  /*
   * print z-values of the geometry and solved density for visualization
   */
  fprintf(g, "CELL_DATA %d\n", nf);
  fprintf(g, "VECTORS normals FLOAT\n");
  for (i = 0; i < nf; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", normals[i].x, normals[i].y,
            normals[i].z);
  }

  fprintf(g, "\n");
  fprintf(g, "VECTORS ds FLOAT\n");
  for (i = 0; i < nf; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", dx[i].x, dx[i].y, dx[i].z);
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS dt FLOAT\n");
  for (i = 0; i < nf; ++i) {
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", dy[i].x, dy[i].y, dy[i].z);
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS Efield_tangential_re FLOAT\n");
  for (i = 0; i < nf; ++i) {
    Efield(vec2, fmid[i], kappa);
    vec = vector3_mul(vec2[0], normals[i]);
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", vec.x, vec.y, vec.z);
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS Efield_tangential_im FLOAT\n");
  for (i = 0; i < nf; ++i) {
    Efield(vec2, fmid[i], kappa);
    vec = vector3_mul(vec2[1], normals[i]);
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", vec.x, vec.y, vec.z);
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS EfieldCurl_tangential_re FLOAT\n");
  for (i = 0; i < nf; ++i) {
    EfieldCurl(vec2, fmid[i], kappa);
    vec = vector3_mul(vec2[0], normals[i]);
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", vec.x, vec.y, vec.z);
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS EfieldCurl_tangential_im FLOAT\n");
  for (i = 0; i < nf; ++i) {
    EfieldCurl(vec2, fmid[i], kappa);
    vec = vector3_mul(vec2[1], normals[i]);
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", vec.x, vec.y, vec.z);
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS Hfield_tangential_re FLOAT\n");
  for (i = 0; i < nf; ++i) {
    Hfield(vec2, fmid[i], kappa);
    vec = vector3_mul(vec2[0], normals[i]);
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", vec.x, vec.y, vec.z);
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS Hfield_tangential_im FLOAT\n");
  for (i = 0; i < nf; ++i) {
    Hfield(vec2, fmid[i], kappa);
    vec = vector3_mul(vec2[1], normals[i]);
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", vec.x, vec.y, vec.z);
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS tested_Re FLOAT\n");
  for (i = 0; i < nf; ++i) {
    vec = vector3_add(vector3_Smul(rhovals[i], dx[i]),
                      vector3_Smul(rhovals[i + 2 * nf], dy[i]));
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", vec.x, vec.y, vec.z);
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS tested_Im FLOAT\n");
  for (i = 0; i < nf; ++i) {
    vec = vector3_add(vector3_Smul(rhovals[i + nf], dx[i]),
                      vector3_Smul(rhovals[i + 3 * nf], dy[i]));
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", vec.x, vec.y, vec.z);
  }
  fprintf(g, "\n");
  fprintf(g, "SCALARS surface_measure FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (i = 0; i < nf; ++i) {
    fprintf(g, "%g\n", surfmeas[i]);
  }
  fprintf(g, "\n");
  fprintf(g, "SCALARS surface_diverence_test_Re FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (i = 0; i < nf; ++i) {
    fprintf(g, "%g\n", (rhodvals[i] + rhodvals[i + 2 * nf]) / surfmeas[i] / h);
  }
  fprintf(g, "\n");
  fprintf(g, "SCALARS surface_diverence_test_Im FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (i = 0; i < nf; ++i) {
    fprintf(g, "%g\n",
            (rhodvals[i + nf] + rhodvals[i + 3 * nf]) / surfmeas[i] / h);
  }
  fprintf(g, "\n");
  fprintf(g, "SCALARS Efield_normal_Re FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (i = 0; i < nf; ++i) {
    Efield(vec2, fmid[i], kappa);
    fprintf(g, "%20.16f\n", vector3_skalp(vec2[0], normals[i]));
  }
  fprintf(g, "\n");
  fprintf(g, "SCALARS Efield_normal_Im FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (i = 0; i < nf; ++i) {
    Efield(vec2, fmid[i], kappa);
    fprintf(g, "%20.16f\n", vector3_skalp(vec2[1], normals[i]));
  }
  fprintf(g, "\n");
  fprintf(g, "SCALARS EfieldCurl_normal_Re FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (i = 0; i < nf; ++i) {
    EfieldCurl(vec2, fmid[i], kappa);
    fprintf(g, "%20.16f\n", vector3_skalp(vec2[0], normals[i]));
  }
  fprintf(g, "\n");
  fprintf(g, "SCALARS EfieldCurl_normal_Im FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (i = 0; i < nf; ++i) {
    EfieldCurl(vec2, fmid[i], kappa);
    fprintf(g, "%20.16f\n", vector3_skalp(vec2[1], normals[i]));
  }
  fprintf(g, "\n");
  fprintf(g, "SCALARS Hfield_normal_Re FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (i = 0; i < nf; ++i) {
    Hfield(vec2, fmid[i], kappa);
    fprintf(g, "%20.16f\n", vector3_skalp(vec2[0], normals[i]));
  }
  fprintf(g, "\n");
  fprintf(g, "SCALARS Hfield_normal_Im FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (i = 0; i < nf; ++i) {
    Hfield(vec2, fmid[i], kappa);
    fprintf(g, "%20.16f\n", vector3_skalp(vec2[1], normals[i]));
  }
  fprintf(g, "\n");
  fprintf(g, "SCALARS patch_ID FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (i = 0; i < nf; ++i) {
    fprintf(g, "%d\n", E[i].patch);
  }
  fprintf(g, "\n");
  fprintf(g, "SCALARS patch_normal_re FLOAT\n");
  fprintf(g, "LOOKUP_TABLE default\n");
  for (i = 0; i < nf; ++i) {
    vec = vector3_add(vector3_Smul(rhovals[i + 0 * nf], dx[i]),
                      vector3_Smul(rhovals[i + 2 * nf], dy[i]));
    if (fabs(E[i].index_s - n / 2) > n / 2 - 4 ||
        fabs(E[i].index_t - n / 2) > n / 2 - 4) {
      if (E[i].index_s <= E[i].index_t && E[i].index_t <= n - E[i].index_s) {
        // s = 0 edge
        vecc = vector3_mul(dy[i], normals[i]);
        vecc = vector3_Smul(1. / vector3_norm(vecc), vecc);
        fprintf(g, "%20.16f\n", vector3_skalp(vec, vecc));
      } else if (E[i].index_s > E[i].index_t &&
                 E[i].index_t <= n - E[i].index_s) {
        // t = 0 edge
        vecc = vector3_mul(normals[i], dx[i]);
        vecc = vector3_Smul(1. / vector3_norm(vecc), vecc);
        fprintf(g, "%20.16f\n", vector3_skalp(vec, vecc));
      } else if (E[i].index_s <= E[i].index_t &&
                 E[i].index_t > n - E[i].index_s) {
        // t = 1 edge
        vecc = vector3_mul(dx[i], normals[i]);
        vecc = vector3_Smul(1. / vector3_norm(vecc), vecc);
        fprintf(g, "%20.16f\n", vector3_skalp(vec, vecc));
      } else if (E[i].index_s > E[i].index_t &&
                 E[i].index_t > n - E[i].index_s) {
        // s = 1 edge
        vecc = vector3_mul(normals[i], dy[i]);
        vecc = vector3_Smul(1. / vector3_norm(vecc), vecc);
        fprintf(g, "%20.16f\n", vector3_skalp(vec, vecc));
      } else {
        fprintf(g, "0\n");
      }
    } else {
      fprintf(g, "0\n");
    }
  }
  fprintf(g, "\n");

  fclose(g);
  free(fmid);
  free(dx);
  free(dy);
  free(normals);
  free(surfmeas);
  free(bigrho);
  return 0;
}

int print_potential_Maxwell(double *Pot, vector3 *Q, int nq, char *dname,
                            double kappa[2]) {
  int i = 0;
  vector3 vec;
  // vector3 f[2];
  vector3 f2[2];
  /*
   * for edge points
   */
  FILE *g = NULL;

  /*
   * visualize the surface of the geometry in vtk format
   */
  g = fopen(dname, "w");

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
  fprintf(g, "POINTS %d FLOAT\n", nq);
  for (i = 0; i < nq; ++i)
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", Q[i].x, Q[i].y, Q[i].z);
  fprintf(g, "\n");

  /*
   * print element list
   */
  fprintf(g, "CELLS %d %d\n", nq, 2 * nq);
  for (i = 0; i < nq; ++i) fprintf(g, "%d %d\n", 1, i);
  fprintf(g, "\n");

  fprintf(g, "CELL_TYPES %d\n", nq);
  for (i = 0; i < nq; ++i) fprintf(g, "%d\n", 1);
  fprintf(g, "\n");

  /*
   * print z-values of the geometry and solved density for visualization
   */
  fprintf(g, "POINT_DATA %d\n", nq);
  fprintf(g, "VECTORS Scattered_me FLOAT\n");
  for (i = 0; i < nq; ++i) {
    vec = vector3_make(Pot[i], Pot[i + 2 * nq], Pot[i + 4 * nq]);
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", vec.x, vec.y, vec.z);
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS Efield_ref FLOAT\n");
  for (i = 0; i < nq; ++i) {
    Efield(f2, Q[i], kappa);
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", f2[0].x, f2[0].y, f2[0].z);
  }
  fprintf(g, "\n");
  fprintf(g, "VECTORS Hfield_ref FLOAT\n");
  for (i = 0; i < nq; ++i) {
    Hfield(f2, Q[i], kappa);
    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", f2[0].x, f2[0].y, f2[0].z);
  }
  fprintf(g, "\n");
  // fprintf(g, "VECTORS Efield_ref FLOAT\n");
  // for (i = 0; i < nq; ++i)
  // {
  //    MaxwellMieIncidentWave(f, Q[i], kappa);
  //    MaxwellMieScatteredWave(f2, Q[i], kappa);
  //    fprintf(g, "%20.16f\t%20.16f\t%20.16f\n", f[0].x + f2[0].x,
  //            f[0].y + f2[0].y, f[0].z + f2[0].z);
  // }
  // fprintf(g, "\n");

  fclose(g);
  return 0;
}

/*-----------------------------------------------------------------------------+
|     int print_geometry                                                       |
|     - visualizes the geometry and the solution in vtk format                 |
|       on the specified level l of the discretization                         |
|     - solution is averaged on the compound panels                            |
|     - might be much slower then the print_geometry routine                   |
|     - prints curved elements!!! (very cool)                                  |
+-----------------------------------------------------------------------------*/
int print_geometry_level(et_node *E, vector3 *P, int nf, double *rho, int l,
                         char *dname) {
  int i = 0;
  int j = 0;
  int max = 0;
  int max_lvl = 0;
  int ppblock = 0;
  int index_s = 0;
  int index_t = 0;
  double dvalue = 0;
  double sc = 1.001;
  FILE *f = NULL;

  /*
   * get pointer to the leafs of the element tree and determine maximum level
   */
  while (E->son[0]) {
    ++max_lvl;
    E = E->son[0];
  }

  /*
   * compute the size of the compound panels
   */
  ppblock = 1 << (max_lvl - l);

  /*
   * get length of pointlist P
   */
  for (i = 0; i < nf; ++i)
    for (j = 0; j < 4; ++j)
      if (E[i].vertex[j] > max) max = E[i].vertex[j];
  ++max;

  /*
   * visualize the surface of the geometry
   */
  f = fopen(dname, "w");

  /*
   * vtk format header
   */
  fprintf(f, "# vtk DataFile Version 3.1\n");
  fprintf(f, "this file hopefully represents my surface now\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  /*
   * print point list
   */
  fprintf(f, "POINTS %d FLOAT\n", max);
  for (i = 0; i < max; ++i)
    fprintf(f, "%20.16f\t%20.16f\t%20.16f\n", sc * P[i].x, sc * P[i].y,
            sc * P[i].z);
  fprintf(f, "\n");

  /*
   * print element list
   */
  fprintf(f, "CELLS %d %d\n", nf / ppblock / ppblock,
          nf / ppblock / ppblock * (ppblock * 4 + 1));

  /*
   * this is the tricky and very slow part: all panels in the compound panel
   * are travelled 4 times to determine the boundary curve of the compound
   * panel
   */
  for (i = 0; i < nf; i += ppblock * ppblock) {
    index_s = E[i].index_s;
    index_t = E[i].index_t;
    fprintf(f, "%d ", 4 * ppblock);
    for (j = 0; j < ppblock * ppblock; ++j)
      if (E[i + j].index_t == index_t) fprintf(f, "%d ", E[i + j].vertex[0]);
    for (j = 0; j < ppblock * ppblock; ++j)
      if (E[i + j].index_s == index_s + ppblock - 1)
        fprintf(f, "%d ", E[i + j].vertex[1]);
    for (j = 0; j < ppblock * ppblock; ++j)
      if (E[i + j].index_t == index_t + ppblock - 1)
        fprintf(f, "%d ", E[i + j].vertex[2]);
    for (j = ppblock * ppblock - 1; j >= 0; --j)
      if (E[i + j].index_s == index_s) fprintf(f, "%d ", E[i + j].vertex[3]);
    fprintf(f, "\n");
  }
  fprintf(f, "\n");

  fprintf(f, "CELL_TYPES %d\n", nf / ppblock / ppblock);
  for (i = 0; i < nf / ppblock / ppblock; ++i) fprintf(f, "%d\n", 7);
  fprintf(f, "\n");

  /*
   * print z-values of the geometry and solved density for visualization
   */
  fprintf(f, "POINT_DATA %d\n", max);
  fprintf(f, "SCALARS Z-value FLOAT\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  for (i = 0; i < max; ++i) fprintf(f, "%20.16f\n", P[i].z);
  fprintf(f, "\n");

  /*
   * unfortunately we have to average the solution
   */
  fprintf(f, "CELL_DATA %d\n", nf / ppblock / ppblock);
  fprintf(f, "SCALARS Cell_Density FLOAT\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  for (i = 0; i < nf; i += ppblock * ppblock) {
    dvalue = 0;
    for (j = 0; j < ppblock * ppblock; ++j) dvalue += rho[i + j];
    dvalue /= ppblock * ppblock;
    fprintf(f, "%20.16f\n", dvalue);
  }

  fclose(f);
  return 0;
}
}  // namespace Bembel