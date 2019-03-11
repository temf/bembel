// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "gauss_square.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "assert.h"
#include "cubature.h"
#include "gauss_legendre.h"
#include "myvector.h"
#include "quadrature.h"

namespace Bembel {
void Tensor_Regel_Square(cubature *Q, quadrature *R)
/*
 * bastelt die Tensor-Produkt-Regel zusammen
 */
// cubature *Q;
// quadrature *R;
{
  int k;

  Q->nop = R->nop * R->nop;
  Q->xi = (vector2 *)malloc(Q->nop * sizeof(vector2));
  Q->w = (double *)malloc(Q->nop * sizeof(double));

  for (k = 0; k < Q->nop; k++) {
    Q->xi[k] = vector2_make(R->xi[k / R->nop], R->xi[k % R->nop]);
    Q->w[k] = R->w[k / R->nop] * R->w[k % R->nop];
  }
  return;
}

void init_Gauss_Square(cubature **Q, int g)
// cubature **Q;
// int g;
{
  quadrature *R;
  int k;

  /*
   * Fehler-Routine
   */

  assert(g < 40 && "gauss_square.c: Quadrature order not supported");
  /*
   * Tensor-Produkte zusammenbasteln
   */
  init_Gauss_Legendre(&R, g);
  (*Q) = (cubature *)malloc(g * sizeof(cubature));
  for (k = 0; k < g; k++) Tensor_Regel_Square(&(*Q)[k], &R[k]);
  free_Gauss_Legendre(&R, g);

  return;
}

void free_Gauss_Square(cubature **Q, int g)
// cubature **Q;
// int g;
{
  int k;

  for (k = 0; k < g; k++) {
    free((*Q)[k].xi);
    free((*Q)[k].w);
  }
  free(*Q);

  return;
}
}  // namespace Bembel