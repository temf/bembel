// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_VECTOR3__
#define __BEMBEL_C_VECTOR3__
/***************
 *  vector3.h  *
 ***************/

/*====================================================*
 *  Kleine Arithmetik fuer dreidimensionale Vektoren  *
 *====================================================*/

#include <math.h>
namespace Bembel {
struct vector2 {
  double x;
  double y;
};

struct vector3 {
  double x;
  double y;
  double z;
};

inline vector2 vector2_make(double x, double y)
/*
 * Typkonvertierung: 2xREAL in vector2
 */
{
  vector2 c;
  c.x = x;
  c.y = y;
  return (c);
}

inline vector2 vector2_add(vector2 a, vector2 b)
/*
 * Vektoraddition
 */
{
  vector2 c;
  c.x = a.x + b.x;
  c.y = a.y + b.y;
  return (c);
}

inline vector2 vector2_sub(vector2 a, vector2 b)
/*
 * Vektorsubtraktion
 */
{
  vector2 c;
  c.x = a.x - b.x;
  c.y = a.y - b.y;
  return (c);
}

inline vector2 vector2_Smul(double a, vector2 b)
/*
 * S-Multiplikation
 */
{
  vector2 c;
  c.x = a * b.x;
  c.y = a * b.y;
  return (c);
}

inline double vector2_skalp(vector2 a, vector2 b)
/*
 * Skalarprodukt
 */
{
  return (a.x * b.x + a.y * b.y);
}

inline double vector2_norm(vector2 a)
/*
 * Euklid-Norm
 */
{
  return (sqrt(a.x * a.x + a.y * a.y));
}

inline vector3 vector3_make(double x, double y, double z)
/*
 * Typkonvertierung: 3xREAL in vector3
 */
{
  vector3 c;
  c.x = x;
  c.y = y;
  c.z = z;
  return (c);
}

inline vector3 vector3_add(vector3 a, vector3 b)
/*
 * Vektoraddition
 */
{
  vector3 c;
  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;
  return (c);
}

inline vector3 vector3_sub(vector3 a, vector3 b)
/*
 * Vektorsubtraktion
 */
{
  vector3 c;
  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
  return (c);
}

inline vector3 vector3_mul(vector3 a, vector3 b)
/*
 * Vektormultiplikation
 */
{
  vector3 c;
  c.x = a.y * b.z - a.z * b.y;
  c.y = a.z * b.x - a.x * b.z;
  c.z = a.x * b.y - a.y * b.x;
  return (c);
}

inline vector3 vector3_Smul(double a, vector3 b)
/*
 * S-Multiplikation
 */
{
  vector3 c;
  c.x = a * b.x;
  c.y = a * b.y;
  c.z = a * b.z;
  return (c);
}

inline double vector3_skalp(vector3 a, vector3 b)
/*
 * Skalarprodukt
 */
{
  return (a.x * b.x + a.y * b.y + a.z * b.z);
}

inline double vector3_norm(vector3 a)
/*
 * Euklid-Norm
 */
{
  return (sqrt(a.x * a.x + a.y * a.y + a.z * a.z));
}

inline vector3 vector3_cross(vector3 a, vector3 b) {
  return vector3_make(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                      a.x * b.y - a.y * b.x);
}
}  // namespace Bembel
#include "mycblas.h"
#endif
