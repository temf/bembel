// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_BEMBDATA_
#define __BEMBEL_C_BEMBDATA_

#include <math.h>
#include <cassert>
#include <complex>
#include "constants.h"
#include "myvector.h"

namespace Bembel {
void spherical(double *z, int m, int n, vector3 x);
double f(vector3 x);
vector3 df(vector3 x);

void Helmholtzf(double d[2], vector3 a, double kappa[2]);
void Efield(vector3 E[2], vector3 x, double kappa[2]);
void EfieldCurl(vector3 E[2], vector3 x, double kappa[2]);
void Efield_minus(vector3 E[2], vector3 x, double kappa[2]);
void Hfield(vector3 H[2], vector3 x, double kappa[2]);
void Hfield_minus(vector3 E[2], vector3 x, double kappa[2]);

void ElectricPlaneWave(vector3 d[2], vector3 a, double kappa[2]);
void ElectricPlaneWaveGrad(vector3 d[2][3], vector3 a, double kappa[2]);
}  // namespace Bembel
#endif
