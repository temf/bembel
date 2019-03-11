// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "gamma.h"
#include "myvector.h"
#include "topology.h"

namespace Bembel {

void init_grid(vector3 **P, int ***F, int m, int *np, int *nf);

void refine_grid(vector3 **P, int ***F, int m, int *np, int *nf);

void init_grid(vector3 **P, int ***F, int m, int *np, int *nf,
               const geometry &Chi)
/*
 * Erstellt die Punkte- und Indexliste fuer Level 0.
 */
// vector3 **P; /* Zeiger auf die Punkteliste */
// int ***F;        /* Zeiger auf die lokale Basisliste */
// int m;           /* 2^m*2^m Patches pro Parametergebiet */
// int *np, *nf;    /* sizeof(P) bzw. sizeof(F) */
{
  const int p = Chi.size();
  int i; /* Laufindizes fuer Parametergebiet */

  /*
   * Speicherplatz allokieren
   */
  (*P) = (vector3 *)malloc(4 * p * sizeof(vector3));
  (*F) = (int **)malloc(p * sizeof(int *));

  /*
   * Eckpunkte der Parametergebiete bestimmen
   */
  *np = 4 * p;
  *nf = p;
  for (i = 0; i < p; i++) {
    /*
     * Bestimme die 4 Eckpunkte des Patches
     */
    (*P)[4 * i] = Chi[i].f(vector2_make(0.0, 0.0));
    (*P)[4 * i + 1] = Chi[i].f(vector2_make(1.0, 0.0));
    (*P)[4 * i + 2] = Chi[i].f(vector2_make(1.0, 1.0));
    (*P)[4 * i + 3] = Chi[i].f(vector2_make(0.0, 1.0));

    (*F)[i] = (int *)malloc(4 * sizeof(int));
    (*F)[i][0] = 4 * i;
    (*F)[i][1] = 4 * i + 1;
    (*F)[i][2] = 4 * i + 2;
    (*F)[i][3] = 4 * i + 3;
  }

  return;
}

void refine_grid(vector3 **P, int ***F, int m, int *np, int *nf,
                 const geometry &Chi)
/*
 * Erstellt die Punkte- und Indexliste fuer alle zusaetzlichen Gitterpunkte
 * des Level m.
 */
// vector3 **P; /* Zeiger auf die Punkteliste */
// int ***F;        /* Zeiger auf die Indexliste */
// int m;           /* 2^m*2^m Patches pro Parametergebiet */
// int *np, *nf;    /* sizeof(P) bzw. sizeof(F) */
{
  const int p = Chi.size();
  int n = 1 << (m - 1); /* n = 2^(m-1) */
  double h = 0.5 / n;   /* Schrittweite auf Level m */
  int i1, i2, i3;       /* Laufindizes */
  int fz;               /* Patchzaehler */

  /*
   * Speicherplatz allokieren
   */
  (*P) = (vector3 *)realloc(*P, 4 * (*np) * sizeof(vector3));
  (*F) = (int **)realloc(*F, 4 * (*nf) * sizeof(int *));
  for (i1 = *nf; i1 < 4 * (*nf); i1++) (*F)[i1] = (int *)calloc(4, sizeof(int));

  /*
   * Kopiere altes Gitter in neues Gitter
   */
  fz = (p - 1) * n * n;
  for (i1 = p - 1; i1 >= 0; i1--) {
    for (i2 = n - 1; i2 >= 0; i2--) {
      for (i3 = n - 1; i3 >= 0; i3--) {
        (*F)[4 * fz + 2 * n * (2 * i2 + 1) + 2 * i3][3] =
            (*F)[fz + n * i2 + i3][3];
        (*F)[4 * fz + 2 * n * (2 * i2 + 1) + 2 * i3 + 1][2] =
            (*F)[fz + n * i2 + i3][2];
        (*F)[4 * fz + 4 * n * i2 + 2 * i3 + 1][1] = (*F)[fz + n * i2 + i3][1];
        (*F)[4 * fz + 4 * n * i2 + 2 * i3][0] = (*F)[fz + n * i2 + i3][0];
      }
    }
    fz -= n * n;
  }

  /*
   * Bestimme Verfeinerung
   */
  fz = 0;
  *nf *= 4;
  for (i1 = 0; i1 < p; i1++) {
    for (i2 = 0; i2 <= 2 * n; i2++) {
      for (i3 = 0; i3 <= 2 * n; i3++) {
        if ((i2 % 2 == 1) || (i3 % 2 == 1)) {
          (*P)[*np] = Chi[i1].f(vector2_make(h * i3, h * i2));
          if ((i2 < 2 * n) && (i3 < 2 * n)) (*F)[fz + 2 * n * i2 + i3][0] = *np;
          if ((i2 < 2 * n) && (0 < i3)) (*F)[fz + 2 * n * i2 + i3 - 1][1] = *np;
          if ((0 < i2) && (0 < i3))
            (*F)[fz + 2 * n * (i2 - 1) + i3 - 1][2] = *np;
          if ((0 < i2) && (i3 < 2 * n))
            (*F)[fz + 2 * n * (i2 - 1) + i3][3] = *np;
          (*np)++;
        }
      }
    }
    fz += 4 * n * n;
  }

  /*
   * Speicherplatz wieder freigeben
   */
  return;
}

int gennet(vector3 **P, int ***F, int M, const geometry &geom)
/*
 * Erstellt die Punkte- und Patchliste in hierarchischer Weise und liefert
 * als Funktionsergebnis die Laenge der Punkteliste, die vom Geschlecht der
 * Oberflaeche abhaengig ist.
 */
// vector3 **P; /* Zeiger auf die Punkteliste */
// int ***F;        /* Zeiger auf die Patchliste */
// int M;           /* 2^m*2^m Patches pro Parametergebiet */
{
  int m;  /* Laufindex fuer das Level */
  int np; /* Laenge von P */
  int nf; /* Laenge von F */

  init_grid(P, F, M, &np, &nf, geom);
  for (m = 1; m <= M; m++) refine_grid(P, F, m, &np, &nf, geom);
  return (np);
}

void unify(vector3 *d, double *r, vector3 d1, double r1, vector3 d2, double r2)
/*
 * bildet die Vereinigung K(d,r) = K(d1,r1) \cup K(d2,r2)
 */
/*vector3 *d,
d1, d2;
double *r, r1, r2;*/
{
  vector3 z;
  double norm;

  z.x = d1.x - d2.x;
  z.y = d1.y - d2.y;
  z.z = d1.z - d2.z;
  norm = sqrt(z.x * z.x + z.y * z.y + z.z * z.z);

  if (norm + r2 <= r1) /* K(d2,r2) \subset K(d1,r1) */
  {
    *d = d1;
    *r = r1;
  } else if (norm + r1 <= r2) /* K(d1,r1) \subset K(d2,r2) */
  {
    *d = d2;
    *r = r2;
  } else /* die Vereinigungsmenge ist keine Kugel */
  {
    (*d).x = 0.5 * (d1.x + d2.x + (r1 - r2) / norm * z.x);
    (*d).y = 0.5 * (d1.y + d2.y + (r1 - r2) / norm * z.y);
    (*d).z = 0.5 * (d1.z + d2.z + (r1 - r2) / norm * z.z);
    *r = 0.5 * (r1 + r2 + norm);
  }
  return;
}

void dist(vector3 *P, int **F, vector3 **D, double **R, int m, const int p)
/*
 * aus der Punkteliste P und der Patchliste F wird fuer jedes Patch i eine
 * Kugel K(D[i],R[i]) gebildet, so dass die Patches bei genuegend feiner
 * Diskretisierung in dieser Kugel enthalten sind.
 */
#if 0
    vector3 *P; /* Punktliste */
int **F;        /* Patchliste */
vector3 **D;    /* zu berechnende Kugelmittelpunkte */
double **R;     /* zu berechnende Kugelradien */
int m;          /* (2^m*2^m) Patches pro Parametergebiet */
#endif
{
  int nf;         /* Anzahl der Patches */
  int k;          /* Laufindex fuer Patch */
  vector3 d1, d2; /* Umkreismittelpunkt von Strecke (0,2) bzw.
                   * (1,3) */
  double r1, r2;  /* Umkreisradius von Strecke (0,2) bzw. (1,3) */

  /*
   * Initialisierung
   */
  nf = p * (1 << m) * (1 << m);
  (*D) = (vector3 *)malloc(nf * sizeof(vector3));
  (*R) = (double *)malloc(nf * sizeof(double));

  for (k = 0; k < nf; k++) {
    unify(&d1, &r1, P[F[k][0]], 0, P[F[k][2]], 0);
    unify(&d2, &r2, P[F[k][1]], 0, P[F[k][3]], 0);
    unify(&(*D)[k], &(*R)[k], d1, r1, d2, r2);
  }

  return;
}

void free_patchlist(int ***F, int nf)
/*
 * gibt den Speicherplatz der (nf,4)-(int)-Patchliste F frei
 */
//     int ***F;
// int nf;
{
  int k;
  for (k = 0; k < nf; k++) free((*F)[k]);
  free(*F);
  return;
}
}  // namespace Bembel