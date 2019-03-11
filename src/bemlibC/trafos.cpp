// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "trafos.h"

namespace Bembel {

/**
 *  \brief         Map from the unit square to the element.
 *
 *  \return        Returns c = a+h*b;
 *
 */
vector2 Kappa(vector2 a, vector2 b, double h)
// vector2 a, b;
// double h;
{
  vector2 c;

  c.x = a.x + h * b.x;
  c.y = a.y + h * b.y;

  return c;
}

/**
 *  \brief         Implements the rotations for the Duffy trick.
 *
 *  \param[in]     b_x            x coordinate of the vector to rotate.
 *  \param[in]     b_y            y coordinate of the vector to rotate.
 *  \param[in]     CASE           Code for the rotation.
 *                                 1: Rotation angle pi/2
 *                                 2: Rotation angle pi
 *                                 3: Rotation angle 3*pi/2
 *                                everything else: identity.
 *
 *  \return        Rotated vector.
 *
 */
vector2 Tau(double b_x, double b_y, int CASE)
// double b_x, b_y;
// int CASE;
{
  switch (CASE) {
    case 1:
      return (vector2_make(1 - b_y, b_x)); /* Drehung um pi/2 */
    case 2:
      return (vector2_make(1 - b_x, 1 - b_y)); /* Drehung um pi */
    case 3:
      return (vector2_make(b_y, 1 - b_x)); /* Drehung um 3*pi/2 */
    default:
      return (vector2_make(b_x, b_y)); /* Identitaet */
  }
}

/**
 *  \brief         Investigates if two patches have similarities.
 *
 *  \param[in]     P              Point list.
 *  \param[in]     F1             Defines the first patch.
 *  \param[in]     F2             Defines the first patch.
 *  \param[out]    ind1           Rotation case for Duffy trick for patch 1.
 *  \param[out]    ind2           Rotation case for Duffy trick for patch 2.
 *
 *  \return        See explanation below.
 *
 *  Investigates if the two patches (F1[0],...,F1[3]) and (F2[0],...,F2[3]) have
 *  a common edge (return 3) or a common point (return 4). In these cases the
 *  rotation variables ind1 and ind2 will be set accordingly. If there are no
 *  similarities return 1.
 *
 */
#ifdef _DK
int compare(vector3 *P, int *F1, int *F2, int *ind1, int *ind2)
// vector3 *P; /* Punkteliste */
// int *F1;
// int *F2; /* Arrays mit den Eckpunkte der gegebenen
// * Patches */
// int *ind1;
// int *ind2; // Indizes der Drehungen fuer die Patches
{
  vector3 d; /* Differenzvektor der Elementeckpunkte */

  /*
   * zuerst auf einen gemeinsamen Punkt untersuchen
   */
  for (*ind1 = 0; *ind1 < 4; (*ind1)++) {
    for (*ind2 = 0; *ind2 < 4; (*ind2)++) {
      d.x = P[F1[*ind1]].x - P[F2[*ind2]].x;
      d.y = P[F1[*ind1]].y - P[F2[*ind2]].y;
      d.z = P[F1[*ind1]].z - P[F2[*ind2]].z;
      if (d.x * d.x + d.y * d.y + d.z * d.z < tol_points) {
        /*
         * gemeinsamer Punkt -> untersuche auf gemneinsame Kante
         */
        d.x = P[F1[(*ind1 + 1) % 4]].x - P[F2[(*ind2 + 3) % 4]].x;
        d.y = P[F1[(*ind1 + 1) % 4]].y - P[F2[(*ind2 + 3) % 4]].y;
        d.z = P[F1[(*ind1 + 1) % 4]].z - P[F2[(*ind2 + 3) % 4]].z;
        if (d.x * d.x + d.y * d.y + d.z * d.z < tol_points) {
          /*
           * normaler Fall: zweiter Punkt bei ind1+1 bzw.
           * ind2-1
           */
          *ind2 = (*ind2 + 3) % 4;
          return (3);
        } else if (*ind1 == 0) {
          d.x = P[F1[3]].x - P[F2[(*ind2 + 1) % 4]].x; /* dies ist der
                                                        * Srfaerfall,
                                                        * denn falls
                                                        * der erste
                                                        * gemeinsame
                                                        * Punkt */
          d.y = P[F1[3]].y - P[F2[(*ind2 + 1) % 4]].y; /* fuer ind1 =
                                                        * 0 aittritt,
                                                        * kann der
                                                        * zweite auch
                                                        * * bei ind3 =
                                                        * 3 liegen! */
          d.z = P[F1[3]].z - P[F2[(*ind2 + 1) % 4]].z;
          if (d.x * d.x + d.y * d.y + d.z * d.z < tol_points) {
            *ind1 = 3;
            return (3);
          }
        }
        return (4);
      }
    }
  }
  return (0);
}
#else
int compare(vector3 *P, int *F1, int *F2, int *ind1, int *ind2)
// vector3 *P;
// int *F1;
// int *F2;
// int *ind1;
// int *ind2;
{
  /*
   * zuerst auf einen gemeinsamen Punkt untersuchen
   */
  for (*ind1 = 0; *ind1 < 4; (*ind1)++) {
    for (*ind2 = 0; *ind2 < 4; (*ind2)++) {
      if (F1[*ind1] == F2[*ind2]) { /* gemeinsamer Punkt -> untersuche auf
                                     * gemneinsame Kante */
        if (F1[3] == F2[(*ind2 + 1) % 4]) { /* dies ist der Sonderfall, denn
                                             * falls der erste */
          *ind1 = 3;  /* gemeinsame Punkt fuer ind1 = 0 auftritt, */
          return (3); /* kann der zweite auch bei ind3 = 3 liegen! */
        } else if (F1[(*ind1 + 1) % 4] ==
                   F2[(*ind2 + 3) % 4]) { /* normaler Fall: zweiter Punkt bei
                                           * ind1+1 bzw. ind2-1 */
          *ind2 = (*ind2 + 3) % 4;
          return (3);
        } else
          return (4);
      }
    }
  }

  return (0);
}
#endif
}  // namespace Bembel