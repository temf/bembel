// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information. #include "data.h"

// namespace Bembel {
// /**
//  *  \brief         Computes a spherical harmonic \f$ z[0]+iz[1] = Y_n^m(x)
//  \f$
//  *                 at the point x.
//  *
//  *  \param[out]    z              Result of the computation
//  *  \param[in]     m              See description
//  *  \param[in]     n              See description
//  *  \param[in]     x              Point on the sphere
//  *
//  */
// void spherical(double *z, int m, int n, vector3 x)
// // double *z;
// // int m;
// // int n;
// // vector3 x;
// {
//   double y1[2];
//   double y2[2];

//   if ((m == 0) && (n == 0)) {
//     z[0] = 0.5 / sqrt(M_PI);
//     z[1] = 0;
//   } else if (m == n) {
//     spherical(y1, m - 1, m - 1, x);
//     z[0] = sqrt((2 * m + 1.0) / (2 * m)) * (x.x * y1[0] - x.y * y1[1]);
//     z[1] = sqrt((2 * m + 1.0) / (2 * m)) * (x.x * y1[1] + x.y * y1[0]);
//   } else if (m + 1 == n) {
//     spherical(z, m, m, x);
//     z[0] *= sqrt(2 * m + 3) * x.z;
//     z[1] *= sqrt(2 * m + 3) * x.z;
//   } else {
//     spherical(y1, m, n - 1, x);
//     spherical(y2, m, n - 2, x);
//     z[0] = sqrt((2 * n + 1.0) / ((n - m) * (n + m))) *
//            (sqrt(2 * n - 1) * x.z * y1[0] -
//             sqrt(((n + m - 1.0) * (n - m - 1.0)) / (2 * n - 3)) * y2[0]);
//     z[1] = sqrt((2 * n + 1.0) / ((n - m) * (n + m))) *
//            (sqrt(2 * n - 1) * x.z * y1[1] -
//             sqrt(((n + m - 1.0) * (n - m - 1.0)) / (2 * n - 3)) * y2[1]);
//   }
//   return;
// }

// /**
//  *  \brief         Provides the function values for the analytical solution
//  and
//  *                 the boundary values.
//  *
//  *  \param[in]     a              Point where to evaluate.
//  *
//  */
// double f(vector3 a)  // vector3 a;
// {
//   /*
//    * return(1.0);
//    */
//   /*
//    * return(a.x*a.x + a.x*a.y + a.z - pow(a.z,5) );
//    */
//   return (4 * a.x * a.x - 3 * a.y * a.y - a.z * a.z);
//   /*
//    * return(1/vector3_norm(vector3_make(2,0,2)));
//    */
//   /*
//    * vector3 b = vector3_make(a.x-1.5,a.y,a.z);
//    * return((b.x+2*b.y+4*b.z)*pow(vector3_norm(b),-3));
//    */
//   /*
//    * double norm, z[2]; norm = vector3_skalp(a, a); if (norm != 0)
//    * spherical(z, 0, 2, vector3_Smul(1 / sqrt(norm), a)); else z[0] = 0;
//    * return (z[0] * norm);
//    */
// }

// /**
//  *  \brief         Provides the derivatives of the analytical solution and
//  *                 the boundary values.
//  *
//  *  \param[in]     a              Point where to evaluate.
//  *
//  */
// vector3 df(vector3 a) {
//   /*
//    * return(vector3_make(0,0,0));
//    */
//   /*
//    * return(vector3_make(2*a.x+a.y,a.x,1-5*pow(a.z,4)));
//    */
//   /*
//    * return(vector3_make(8*a.x,-6*a.y,-2*a.z));
//    */
//   /*
//    * vector3 r = vector3_make(2,0,2);
//    * return(vector3_Smul(1/pow(vector3_norm(r),3),r));
//    */
//   /*
//    * vector3 b = vector3_make(a.x,a.y,a.z); double c = vector3_norm(b);
//    double
//    * d = b.x+2*b.y+4*b.z; double e = pow(c,-5);
//    *
//    return(vector3_make((c*c-3*b.x*d)*e,(2*c*c-3*b.y*d)*e,(4*c*c-3*b.z*d)*e));
//    */
//   double norm, z[2];
//   vector3 x;
//   norm = vector3_norm(a);
//   x = vector3_Smul(1 / norm, a);
//   if (norm != 0)
//     spherical(z, 0, 2, x);
//   else
//     z[0] = 0;
//   return (vector3_Smul(z[0] * 2 * norm, x));
// }

// /**
//  *  \brief         Provides the function values for the analytical solution
//  and
//  *                 the boundary values.
//  *
//  *  \param[out]    d              Function output.
//  *  \param[in]     a              Point where to evaluate.
//  *
//  */
// void Helmholtzf(double d[2], vector3 a, double kappa[2]) {
//   double r, e;
//   r = vector3_norm(vector3_make(a.x - 0.1, a.y - 0., a.z - 0.));
//   e = exp(-kappa[1] * r) / r;
//   d[0] = cos(kappa[0] * r) * e;
//   d[1] = sin(kappa[0] * r) * e;
//   return;
// }

// static const vector3 y = vector3_make(0.2, 0.2, 0.2);  // position
// static const vector3 p = vector3_make(0., 0.1, 0.1);   // länge

// /**
//  * \brief         This function implements the electric field of a Hertz
//  dipole.
//  *
//  * \param[out]    E           Electric field.
//  * \param[in]     x           Point where to evaluate.
//  * \param[in]     kappa       Wave number.
//  *
//  * \attention     If the sign in the fundamental solution changes, the sign
//  *                must be changed here, too.
//  *
//  * \author        Juergen Doelz
//  */
// void Efield(vector3 E[2], vector3 x, double kappa[2]) {
//   double r, r2, r3;
//   std::complex<double> imag(0, 1);
//   std::complex<double> expc;
//   std::complex<double> ec;
//   std::complex<double> kappac(-kappa[0], -kappa[1]);
//   std::complex<double> kappac2 = kappac * kappac;
//   vector3 c, h;

//   c.x = x.x - y.x;
//   c.y = x.y - y.y;
//   c.z = x.z - y.z;
//   r = sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
//   r2 = r * r;
//   r3 = r2 * r;
//   vector3 n = vector3_Smul(1. / r, c);

//   h = vector3_mul(vector3_mul(n, p), n);

//   expc = std::exp(imag * kappac * r);

//   /* E =  exp(ikr) * kappa^2*(n x p) x n / r */
//   ec = kappac2 * expc / r;
//   E[0] = vector3_Smul(std::real(ec), h);
//   E[1] = vector3_Smul(std::imag(ec), h);

//   /* E += (3*(n.p)*n-p)*(1/r3*e+i*kappa/r2) */
//   h = vector3_add(vector3_Smul(3 * vector3_skalp(n, p), n),
//                   vector3_Smul(-1., p));
//   ec = (1 / r3 - imag * kappac / r2) * expc;
//   E[0] = vector3_add(E[0], vector3_Smul(std::real(ec), h));
//   E[1] = vector3_add(E[1], vector3_Smul(std::imag(ec), h));
//   return;
// }

// /**
//  * \brief         This function implements the curl of Efield().
//  *
//  * \param[out]    H           Magnetic field.
//  * \param[in]     x           Point where to evaluate.
//  * \param[in]     kappa       Wave number.
//  *
//  * \attention     If the sign in the fundamental solution changes, the sign
//  *                must be changed here, too.
//  *
//  */
// void EfieldCurl(vector3 H[2], vector3 x, double kappa[2]) {
//   double r, r2;
//   std::complex<double> imag(0, 1);
//   std::complex<double> ec;
//   std::complex<double> kappac(-kappa[0], -kappa[1]);
//   std::complex<double> kappac2 = kappac * kappac;
//   vector3 c;

//   c.x = x.x - y.x;
//   c.y = x.y - y.y;
//   c.z = x.z - y.z;
//   r = sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
//   r2 = r * r;
//   vector3 n = vector3_Smul(1. / r, c);

//   /* exp(ikr) */
//   ec = std::exp(imag * kappac * r) * (kappac2 / r - 1. / (imag * kappac *
//   r2));

//   ec *= imag * kappac;

//   c = vector3_mul(n, p);
//   H[0] = vector3_Smul(std::real(ec), c);
//   H[1] = vector3_Smul(std::imag(ec), c);

//   return;
// }

// void Efield_minus(vector3 E[2], vector3 x, double kappa[2]) {
//   Efield(E, x, kappa);
//   E[0] = vector3_Smul(-1., E[0]);
//   E[1] = vector3_Smul(-1., E[1]);
//   return;
// }

// /**
//  * \brief         This function implements the magnetic field of a Hertz
//  dipole.
//  *                It is the curl of the electric field Efield().
//  *
//  * \param[out]    H           Magnetic field.
//  * \param[in]     x           Point where to evaluate.
//  * \param[in]     kappa       Wave number.
//  *
//  * \attention     If the sign in the fundamental solution changes, the sign
//  *                must be changed here, too.
//  *
//  */
// void Hfield(vector3 H[2], vector3 x, double kappa[2]) {
//   double r, r2;
//   std::complex<double> imag(0, 1);
//   std::complex<double> ec;
//   std::complex<double> kappac(kappa[0], kappa[1]);
//   std::complex<double> kappac2 = kappac * kappac;
//   vector3 c;

//   c.x = x.x - y.x;
//   c.y = x.y - y.y;
//   c.z = x.z - y.z;
//   r = sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
//   r2 = r * r;
//   vector3 n = vector3_Smul(1. / r, c);

//   /* exp(ikr) */
//   ec = std::exp(imag * kappac * r) * (kappac2 / r - 1. / (imag * kappac *
//   r2));

//   c = vector3_mul(n, p);
//   H[0] = vector3_Smul(std::real(ec), c);
//   H[1] = vector3_Smul(std::imag(ec), c);

//   return;
// }

// void Hfield_minus(vector3 E[2], vector3 x, double kappa[2]) {
//   Hfield(E, x, kappa);
//   E[0] = vector3_Smul(-1., E[0]);
//   E[1] = vector3_Smul(-1., E[1]);
//   return;
// }

// void ElectricPlaneWave(vector3 d[2], vector3 a, double kappa[2]) {
//   assert(!kappa[1]);
//   d[0].x = cos(-kappa[0] * a.z);
//   d[0].y = 0;
//   d[0].z = 0;
//   d[1].x = sin(-kappa[0] * a.z);
//   d[1].y = 0;
//   d[1].z = 0;
//   return;
// }

// void ElectricPlaneWaveCurl(vector3 d[2], vector3 a, double kappa[2]) {
//   assert(!kappa[1]);
//   d[0].x = 0.;
//   d[0].y = kappa[0] * sin(-kappa[0] * a.z);
//   d[0].z = 0.;
//   d[1].x = 0.;
//   d[1].y = -kappa[0] * cos(-kappa[0] * a.z);
//   d[1].z = 0.;
//   return;
// }

// void ElectricPlaneWaveGrad(vector3 d[2][3], vector3 a, double kappa[2]) {
//   assert(!kappa[1]);
//   d[0][0].x = 0;
//   d[0][0].y = 0;
//   d[0][0].z = kappa[0] * sin(-kappa[0] * a.z);
//   d[0][1].x = 0;
//   d[0][1].y = 0;
//   d[0][1].z = 0;
//   d[0][2].x = 0;
//   d[0][2].y = 0;
//   d[0][2].z = 0;
//   d[1][0].x = 0;
//   d[1][0].y = 0;
//   d[1][0].z = -kappa[0] * cos(-kappa[0] * a.z);
//   d[1][1].x = 0;
//   d[1][1].y = 0;
//   d[1][1].z = 0;
//   d[1][2].x = 0;
//   d[1][2].y = 0;
//   d[1][2].z = 0;
//   return;
// }

// /**
//  * \brief Gradient of Curl
//  */
// void ElectricPlaneWaveGradCurl(vector3 d[2][3], vector3 a, double kappa[2]) {
//   assert(!kappa[1]);
//   d[0][0].x = 0;
//   d[0][0].y = 0;
//   d[0][0].z = 0;
//   d[0][1].x = 0;
//   d[0][1].y = 0;
//   d[0][1].z = -kappa[0] * kappa[0] * cos(-kappa[0] * a.z);
//   d[0][2].x = 0;
//   d[0][2].y = 0;
//   d[0][2].z = 0;
//   d[1][0].x = 0;
//   d[1][0].y = 0;
//   d[1][0].z = 0;
//   d[1][1].x = 0;
//   d[1][1].y = 0;
//   d[1][1].z = -kappa[0] * kappa[0] * sin(-kappa[0] * a.z);
//   d[1][2].x = 0;
//   d[1][2].y = 0;
//   d[1][2].z = 0;
//   return;
//   return;
// }

// void MieSurfaceCurrentIndirect(vector3 d[2], vector3 a, double kappa[2]) {
// #ifndef _USE_MIE_
//   assert(!"Compile with -D_USE_MIE_ to enable Mie scattering reference
//   solution\n");
// #else
//   assert(!kappa[1]);
//   mie_current_indirect(kappa[0], a, d);
// #endif
//   return;
// }

// void MieScatteredWave(vector3 d[2], vector3 a, double kappa[2]) {
// #ifndef _USE_MIE_
//   assert(!"Compile with -D_USE_MIE_ to enable Mie scattering reference
//   solution\n");
// #else
//   assert(!kappa[1]);
//   mie_field(kappa[0], a, d);
// #endif
//   return;
// }

// void MieScatteredWaveRadius(vector3 d[2], vector3 a, double kappa[2],
//                             double r) {
//   double kappar[2] = {kappa[0] * r, kappa[1] * r};
//   vector3 ar = vector3_Smul(1 / r, a);
//   MieScatteredWave(d, ar, kappar);
//   return;
// }
// }  // namespace Bembel