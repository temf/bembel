// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "init_gridpoints.h"

namespace Bembel {
/**
 * \brief Generates a grid of cubes containing the loaded geometry where outer
 * cubes are marked by -1 and cubes containing the boundary are marked by 1. All
 * other cubes are marked by 0.
 *
 * \param[out]   cubes:   Array of cubes of length size_x*size_y*size_z
 * \param[in]    h:       Edge length of the cubes
 * \param[out]   size_x, size_y, size_z:
 *                        Edge length of the grid in direction x,y,z
 * \param[out]   xmin, ymin, zmin, xmax, ymax, zmax:
 *                        Minimum/maximum coordinate of the grid edges in
 * direction x,y,z
 *
 * The algorithm implements a heuristic to generate sufficient enough points on
 * the surface to be almost sure that the distance of the generated points in
 * space is much less than h. Those points are then used to determine the
 * minimum and the maximum coordinate of the surface in x,y,z direction and to
 * generate a grid containing the whole geometry. Outer cubes of the grid are
 * marked by -1, cubes containing the boundary of the geometry are marked by 1.
 * All other cubes are marked by 0.
 *
 * \warning The implemented heuristic computes enough points to work well with
 * the usual geometries. However it may fail if one patch has a quickly changing
 *          parametrization.
 *
 */
int get_cube_grid(int **cubes, double h, int *size_x, int *size_y, int *size_z,
                  double *xmin, double *ymin, double *zmin, double *xmax,
                  double *ymax, double *zmax, const geometry &Chi) {
  int i1, i2, i3;
  const int p = Chi.size();
  // int p = init_p(); /* number of patches */
  int np = 0; /* number of points in one direction to
               * generate on the current patch */
  double xmean = 0, ymean = 0, zmean = 0;
  double mel = 0;
  vector3 corner1, corner2, corner3, corner4;
  vector3 cornerdiff1, cornerdiff2, cornerdiff3, cornerdiff4;
  vector3 *Plist;
  vector3 *Nlist;
  int xcoord, ycoord, zcoord;
  int Plist_size = 0; /* size of Plist and Nlist */
  // parametrix *Chi = NULL;

  // init_Chi(&Chi);

  Plist = (vector3 *)calloc(1, sizeof(vector3));
  Nlist = (vector3 *)calloc(1, sizeof(vector3));

  for (i1 = 0; i1 < p; ++i1) {
    /*
     * compute coordinates of corners
     */
    corner1 = Chi[i1].f(vector2_make(0.0, 0.0));
    corner2 = Chi[i1].f(vector2_make(1.0, 0.0));
    corner3 = Chi[i1].f(vector2_make(1.0, 1.0));
    corner4 = Chi[i1].f(vector2_make(0.0, 1.0));

    /*
     * compute differences between corners
     */
    cornerdiff1 = vector3_sub(corner1, corner2);
    cornerdiff2 = vector3_sub(corner2, corner3);
    cornerdiff3 = vector3_sub(corner3, corner4);
    cornerdiff4 = vector3_sub(corner4, corner1);

    /*
     * determine maximal distance of corners
     */
    mel = 0;
    mel = std::max(mel, vector3_norm(cornerdiff1));
    mel = std::max(mel, vector3_norm(cornerdiff2));
    mel = std::max(mel, vector3_norm(cornerdiff3));
    mel = std::max(mel, vector3_norm(cornerdiff4));

    /*
     * determine number of points in each direction
     */
    np = 16 * (int)(mel / h);

    Plist = (vector3 *)realloc(Plist, (Plist_size + np * np) * sizeof(vector3));
    Nlist = (vector3 *)realloc(Nlist, (Plist_size + np * np) * sizeof(vector3));

    for (i2 = 0; i2 < np; ++i2)
      for (i3 = 0; i3 < np; ++i3) {
        Plist[Plist_size + np * i2 + i3] =
            Chi[i1].f(vector2_make((double)(i2 + 1) / (double)(np + 1),
                                   (double)(i3 + 1) / (double)(np + 1)));
        Nlist[Plist_size + np * i2 + i3] =
            Chi[i1].n_f(vector2_make((double)(i2 + 1) / (double)(np + 1),
                                     (double)(i3 + 1) / (double)(np + 1)));
      }

    Plist_size += np * np;
  }

  *xmin = Plist[0].x;
  *xmax = Plist[0].x;
  *ymin = Plist[0].y;
  *ymax = Plist[0].y;
  *zmin = Plist[0].z;
  *zmax = Plist[0].z;

  for (i1 = 1; i1 < Plist_size; ++i1) {
    *xmin = std::min(*xmin, Plist[i1].x);
    *xmax = std::max(*xmax, Plist[i1].x);
    *ymin = std::min(*ymin, Plist[i1].y);
    *ymax = std::max(*ymax, Plist[i1].y);
    *zmin = std::min(*zmin, Plist[i1].z);
    *zmax = std::max(*zmax, Plist[i1].z);
  }

  xmean = (*xmin + *xmax) / 2;
  ymean = (*ymin + *ymax) / 2;
  zmean = (*zmin + *zmax) / 2;

  *size_x = (int)ceil((*xmax - xmean) / h) + 1;
  *size_y = (int)ceil((*ymax - ymean) / h) + 1;
  *size_z = (int)ceil((*zmax - zmean) / h) + 1;

  *xmin = xmean - h * (*size_x);
  *ymin = ymean - h * (*size_y);
  *zmin = zmean - h * (*size_z);

  *xmax = xmean + h * (*size_x) + h;
  *ymax = ymean + h * (*size_y) + h;
  *zmax = zmean + h * (*size_z) + h;

  (*size_x) *= 2;
  (*size_y) *= 2;
  (*size_z) *= 2;

  (*size_x) += 1;
  (*size_y) += 1;
  (*size_z) += 1;

  *cubes = (int *)calloc((*size_x) * (*size_y) * (*size_z), sizeof(int));

  /*************************************************************************************/
  /*
   * mark boundary cubes as outside with -1
   */

  /*************************************************************************************/
  for (i1 = 0; i1 < (*size_x) * (*size_y); ++i1) {
    (*cubes)[i1] = -1;
    (*cubes)[(*size_x) * (*size_y) * (*size_z) - i1 - 1] = -1;
  }
  for (i1 = 0; i1 < (*size_x); ++i1) {
    for (i2 = 0; i2 < (*size_z); ++i2) {
      (*cubes)[i1 + (*size_x) * (*size_y) * i2] = -1;
      (*cubes)[(*size_x) * (*size_y) * (*size_z) - i1 -
               (*size_x) * (*size_y) * i2 - 1] = -1;
    }
  }
  for (i1 = 0; i1 < (*size_y); ++i1) {
    for (i2 = 0; i2 < (*size_z); ++i2) {
      (*cubes)[(*size_x) * i1 + (*size_x) * (*size_y) * i2] = -1;
      (*cubes)[(*size_x) * (*size_y) * (*size_z) - (*size_x) * i1 -
               (*size_x) * (*size_y) * i2 - 1] = -1;
    }
  }

  /*************************************************************************************/
  /*
   * mark cubes containing the boundary of the computational domain with 1
   */

  /*************************************************************************************/
  for (i1 = 0; i1 < Plist_size; ++i1) {
    if (Plist[i1].x > xmean) {
      xcoord = (*size_x) - (int)floor(((*xmax) - Plist[i1].x) / h) - 1;
      if (fabs(floor(((*xmax) - Plist[i1].x) / h) -
               ((*xmax) - Plist[i1].x) / h) < 1e-12 &&
          Nlist[i1].x > 0)
        xcoord++;
    } else {
      xcoord = (int)floor((Plist[i1].x - (*xmin)) / h);
      if (fabs(floor((Plist[i1].x - (*xmin)) / h) -
               (Plist[i1].x - (*xmin)) / h) < 1e-12 &&
          Nlist[i1].x < 0)
        xcoord--;
    }

    if (Plist[i1].y > ymean) {
      ycoord = (*size_y) - (int)floor(((*ymax) - Plist[i1].y) / h) - 1;
      if (fabs(floor(((*ymax) - Plist[i1].y) / h) -
               ((*ymax) - Plist[i1].y) / h) < 1e-12 &&
          Nlist[i1].y > 0)
        ycoord++;
    } else {
      ycoord = (int)floor((Plist[i1].y - (*ymin)) / h);
      if (fabs(floor((Plist[i1].y - (*ymin)) / h) -
               (Plist[i1].y - (*ymin)) / h) < 1e-12 &&
          Nlist[i1].y < 0)
        ycoord--;
    }

    if (Plist[i1].z > zmean) {
      zcoord = (*size_z) - (int)floor(((*zmax) - Plist[i1].z) / h) - 1;
      if (fabs(floor(((*zmax) - Plist[i1].z) / h) -
               ((*zmax) - Plist[i1].z) / h) < 1e-12 &&
          Nlist[i1].z > 0)
        zcoord++;
    } else {
      zcoord = (int)floor((Plist[i1].z - (*zmin)) / h);
      if (fabs(floor((Plist[i1].z - (*zmin)) / h) -
               (Plist[i1].z - (*zmin)) / h) < 1e-12 &&
          Nlist[i1].z < 0)
        zcoord--;
    }

    (*cubes)[xcoord + (*size_x) * ycoord + (*size_x) * (*size_y) * zcoord] = 1;
  }

  free(Plist);
  free(Nlist);
  // free_Chi(&Chi);

  return 0;
}

/**
 * \brief Determines if the point P is already contained in the point list Q, if
 * not it will be attached to the end of it. Returns the index of the point in
 *        the point list.
 *
 * \param[in, out]   Q:   Point list
 * \param[in, out]   nq:  Lenth of Q
 * \param[in]        P:   Point which will be searched/attached
 * \return           Index of P in Q
 *
 */
int add_gridpoint(vector3 **Q, int *nq, vector3 P)
// vector3 **Q;
// int *nq;
// vector3 P;
{
  int i;
  for (i = 0; i < *nq; ++i) {
    if (fabs((*Q)[i].x - P.x) < 1e-12 && fabs((*Q)[i].y - P.y) < 1e-12 &&
        fabs((*Q)[i].z - P.z) < 1e-12)
      break;
  }

  if (i == (*nq)) {
    *Q = (vector3 *)realloc(*Q, ((*nq) + 1) * sizeof(vector3));
    (*Q)[*nq] = P;
    (*nq)++;
  }

  return i;
}

/**
 * \brief Fills the loaded geometry with a grid of cubes with edge length h.
 *
 * \param[out]   Q:      Array in which the edges of the cubes are stored
 * \param[out]   nq:     Length of Q
 * \param[out]   Cells:  Array which describes which indices in Q form a Cube.
 * The points in Q are taken as edges. \param[out]   nc:     Length of Cells
 * \param[in]    h:      Length of the cube edges
 * \param[in]    h_disc: Minimum distance of the cube to the boundary
 *
 * This function takes a grid which is generated by get_cube_grid() and where
 * the outer cells are marked by -1 and the cells containing the boundary by 1.
 * All other cells must be marked by zero. The algorithm then marks all cells
 * marked as 0 and touching a cell marked by -1 by -1 until no changes can be
 * made anymore. All cells still marked by 0 are contained in the loaded
 * geometry. Following this part of the algorithm, the algorithm continues to
 * remove the outer cells still marked by 0 until the distance given by h_disc
 * to the boundary is reached.
 *
 * \warning It is recommended to set h<=h_disc
 *
 * On Harbrecht's geometries good values to generate something between 20'000
 * and 50'000 points are:
 *  * crankshaft: h,h_disc = 0.1
 *  * sphere, toy : h,h_disc = 0.05
 *
 */
int init_gridpoints(vector3 **Q, int *nq, gridcell **Cells, int *nc, double h,
                    double h_disc, const geometry &Chi)
// vector3 **Q;
// int *nq;
// gridcell **Cells;
// int *nc;
// double h;
// double h_disc;
{
  int i1, i2, i3;
  int *cubes = NULL;
  double xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0;
  int size_x, size_y, size_z;
  int bbr = 0; /* number of boundary blocks removed */
  int madeachange = 1;

  /*************************************************************************************/
  /*
   * get boundary points on the grid and determine grid size
   */

  /*************************************************************************************/

  get_cube_grid(&cubes, h, &size_x, &size_y, &size_z, &xmin, &ymin, &zmin,
                &xmax, &ymax, &zmax, Chi);

  /*************************************************************************************/
  /*
   * mark not yet marked cubes touching a cube marked as outside with -1
   */

  /*************************************************************************************/
  while (madeachange == 1) {
    madeachange = 0;
    for (i1 = 1; i1 < size_x - 1; ++i1) {
      for (i2 = 1; i2 < size_y - 1; ++i2) {
        for (i3 = 1; i3 < size_z - 1; ++i3) {
          if (cubes[i1 + size_x * i2 + size_x * size_y * i3] == 0)
            if (cubes[i1 - 1 + size_x * i2 + size_x * size_y * i3] == -1 ||
                cubes[i1 + 1 + size_x * i2 + size_x * size_y * i3] == -1 ||
                cubes[i1 + size_x * (i2 - 1) + size_x * size_y * i3] == -1 ||
                cubes[i1 + size_x * (i2 + 1) + size_x * size_y * i3] == -1 ||
                cubes[i1 + size_x * i2 + size_x * size_y * (i3 - 1)] == -1 ||
                cubes[i1 + size_x * i2 + size_x * size_y * (i3 + 1)] == -1) {
              cubes[i1 + size_x * i2 + size_x * size_y * i3] = -1;
              madeachange = 1;
            }
        }
      }
    }
  }

  /*************************************************************************************/
  /*
   * mark outer cubes of the geometry with incrementing numbers until distance
   * to
   */
  /*
   * boundary is sufficient
   */

  /*************************************************************************************/
  for (bbr = 0; bbr * h < h_disc * sqrt(3); ++bbr) {
    for (i1 = 1; i1 < size_x - 1; ++i1) {
      for (i2 = 1; i2 < size_y - 1; ++i2) {
        for (i3 = 1; i3 < size_z - 1; ++i3) {
          if (cubes[i1 + size_x * i2 + size_x * size_y * i3] == 0 &&
              (cubes[i1 - 1 + size_x * i2 + size_x * size_y * i3] == bbr + 1 ||
               cubes[i1 + 1 + size_x * i2 + size_x * size_y * i3] == bbr + 1 ||
               cubes[i1 + size_x * (i2 - 1) + size_x * size_y * i3] ==
                   bbr + 1 ||
               cubes[i1 + size_x * (i2 + 1) + size_x * size_y * i3] ==
                   bbr + 1 ||
               cubes[i1 + size_x * i2 + size_x * size_y * (i3 - 1)] ==
                   bbr + 1 ||
               cubes[i1 + size_x * i2 + size_x * size_y * (i3 + 1)] == bbr + 1))
            cubes[i1 + size_x * i2 + size_x * size_y * i3] = bbr + 2;
        }
      }
    }
  }

  /*************************************************************************************/
  /*
   * insert points from a cube being marked by 0 into point list
   */

  /*************************************************************************************/
  *nq = 0;
  *Q = (vector3 *)calloc(1, sizeof(vector3));
  *nc = 0;
  *Cells = (gridcell *)calloc(1, sizeof(gridcell));

  for (i1 = 0; i1 < size_x - 0; ++i1) {
    for (i2 = 0; i2 < size_y - 0; ++i2) {
      for (i3 = 0; i3 < size_z - 0; ++i3) {
        if (cubes[i1 + size_x * i2 + size_x * size_y * i3] == 0) {
          *Cells = (gridcell *)realloc(*Cells, (*nc + 1) * sizeof(gridcell));

          (*Cells)[*nc].c1 = add_gridpoint(
              Q, nq, vector3_make(i1 * h + xmin, i2 * h + ymin, i3 * h + zmin));
          (*Cells)[*nc].c2 = add_gridpoint(
              Q, nq,
              vector3_make((i1 + 1) * h + xmin, i2 * h + ymin, i3 * h + zmin));
          (*Cells)[*nc].c3 =
              add_gridpoint(Q, nq,
                            vector3_make((i1 + 1) * h + xmin,
                                         (i2 + 1) * h + ymin, i3 * h + zmin));
          (*Cells)[*nc].c4 = add_gridpoint(
              Q, nq,
              vector3_make(i1 * h + xmin, (i2 + 1) * h + ymin, i3 * h + zmin));
          (*Cells)[*nc].c5 = add_gridpoint(
              Q, nq,
              vector3_make(i1 * h + xmin, i2 * h + ymin, (i3 + 1) * h + zmin));
          (*Cells)[*nc].c6 =
              add_gridpoint(Q, nq,
                            vector3_make((i1 + 1) * h + xmin, i2 * h + ymin,
                                         (i3 + 1) * h + zmin));
          (*Cells)[*nc].c7 = add_gridpoint(
              Q, nq,
              vector3_make((i1 + 1) * h + xmin, (i2 + 1) * h + ymin,
                           (i3 + 1) * h + zmin));
          (*Cells)[*nc].c8 =
              add_gridpoint(Q, nq,
                            vector3_make(i1 * h + xmin, (i2 + 1) * h + ymin,
                                         (i3 + 1) * h + zmin));

          (*nc)++;
        }
      }
    }
  }

  printf("                           %d points\n", *nq);
  printf("                           %d cells\n", *nc);

  free(cubes);

  return 0;
}

/**
 * \brief Writes a vtk file with points generated by init_gridpoints() and
 * values on this point.
 *
 * \param[in]  P:       Point list, generated by init_gridpoints()
 * \param[in]  nq:      Length of P, generated by init_gridpoints()
 * \param[in]  Pot:     Value of the points contained in P
 * \param[in]  Cells:   Describes which points form a cube together,
 *                      generated by init_gridpoints()
 * \param[in]  nc:      Length of Cells, generated by init_gridpoints()
 * \param[in]  dname:   filename of the output file
 *
 */
int print_cells(vector3 *P, int nq, double *Pot, gridcell *Cells, int nc,
                char *dname) {
  int i = 0;
  FILE *f = NULL;

  /*
   * visualize the surface of the geometry in vtk format
   */
  f = fopen(dname, "w");

  /*
   * vtk format header
   */
  fprintf(f, "# vtk DataFile Version 3.1\n");
  fprintf(f, "this file hopefully represents my potential now\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  /*
   * print point list
   */
  fprintf(f, "POINTS %d FLOAT\n", nq);
  for (i = 0; i < nq; ++i)
    fprintf(f, "%20.16f\t%20.16f\t%20.16f\n", P[i].x, P[i].y, P[i].z);
  fprintf(f, "\n");

  /*
   * print element list
   */
  fprintf(f, "CELLS %d %d\n", nc, 9 * nc);
  for (i = 0; i < nc; ++i)
    fprintf(f, "%d %d %d %d %d %d %d %d %d\n", 8, Cells[i].c1, Cells[i].c2,
            Cells[i].c3, Cells[i].c4, Cells[i].c5, Cells[i].c6, Cells[i].c7,
            Cells[i].c8);
  fprintf(f, "\n");

  fprintf(f, "CELL_TYPES %d\n", nc);
  for (i = 0; i < nc; ++i) fprintf(f, "%d\n", 12);
  fprintf(f, "\n");

  /*
   * print z-values of the geometry and solved density for visualization
   */
  fprintf(f, "POINT_DATA %d\n", nq);
  fprintf(f, "SCALARS Potential FLOAT\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  for (i = 0; i < nq; ++i) fprintf(f, "%20.16f\n", Pot[i]);
  fprintf(f, "\n");

  fprintf(f, "CELL_DATA %d\n", nc);
  fprintf(f, "SCALARS Cell_mean FLOAT\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  for (i = 0; i < nc; ++i)
    fprintf(f, "%20.16f\n",
            (Pot[Cells[i].c1] + Pot[Cells[i].c2] + Pot[Cells[i].c3] +
             Pot[Cells[i].c4] + Pot[Cells[i].c5] + Pot[Cells[i].c6] +
             Pot[Cells[i].c7] + Pot[Cells[i].c8]) /
                8);
  fprintf(f, "\n");

  fprintf(f, "SCALARS Cell_min FLOAT\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  for (i = 0; i < nc; ++i)
    fprintf(f, "%20.16f\n",
            std::min(
                Pot[Cells[i].c1],
                std::min(
                    Pot[Cells[i].c2],
                    std::min(
                        Pot[Cells[i].c3],
                        std::min(
                            Pot[Cells[i].c4],
                            std::min(Pot[Cells[i].c5],
                                     std::min(Pot[Cells[i].c6],
                                              std::min(Pot[Cells[i].c7],
                                                       Pot[Cells[i].c8]))))))));
  fprintf(f, "\n");

  fprintf(f, "SCALARS Cell_max FLOAT\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  for (i = 0; i < nc; ++i)
    fprintf(f, "%20.16f\n",
            std::max(
                Pot[Cells[i].c1],
                std::max(
                    Pot[Cells[i].c2],
                    std::max(
                        Pot[Cells[i].c3],
                        std::max(
                            Pot[Cells[i].c4],
                            std::max(Pot[Cells[i].c5],
                                     std::max(Pot[Cells[i].c6],
                                              std::max(Pot[Cells[i].c7],
                                                       Pot[Cells[i].c8]))))))));

  fclose(f);
  return 0;
}
}  // namespace Bembel