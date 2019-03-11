// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_C_PHI__
#define __BEMBEL_C_PHI__

#include <assert.h>
#include "myvector.h"
#include "spline/bernstein.h"
namespace Bembel {
void phi_0P(double *, double, double);
void phiphi_0P(double *, vector2);
void phiphi_dx_0P(double *, vector2);
void phiphi_dy_0P(double *, vector2);
void Phi_times_Phi_0P(double *, double, vector2, vector2);

void phi_1P(double *, double, double);
void phiphi_1P(double *, vector2);
void phiphi_dx_1P(double *, vector2);
void phiphi_dy_1P(double *, vector2);
void Phi_times_Phi_1P(double *, double, vector2, vector2);
void VPhi_scal_VPhi_1P(double *, double, vector2, vector2, vector3, vector3,
                       vector3, vector3);
void Div_Phi_times_Div_Phi_1P(double *, double, vector2, vector2);
void Curl_Phi_times_Curl_Phi_1P(double *, double, vector2, vector2, vector3,
                                vector3, vector3, vector3);

void phi_2B(double *, double, double);
void phiphi_2B(double *, vector2);
void phiphi_dx_2B(double *, vector2);
void phiphi_dy_2B(double *, vector2);
void Phi_times_Phi_2B(double *, double, vector2, vector2);

// Archetypes, supposed to fail.

// template <int N> inline void phiBase(double *array, const double x) noexcept
// {
//   // assert(N && 0 && "Something Went Wrong. This should not be called -->
//   // phi.h");
//   std::cout << "Phi.h: This is horrible. This should never happen!\n";
//   return;
// }

// template <int N>
// inline void phiBaseDer(double *array, const double x) noexcept {
//   // assert(N && 0 && "Something Went Wrong. This should not be called -->
//   // phi.h");
//   std::cout << "Phi.h: This is horrible. This should never happen!\n";
//   return;
// }
// // Specializations

// //************************************************************************
// //********* PHI
// //************************************************************************

// template <> inline void phiBase<1>(double *array, const double x) noexcept {
//   array[0] = 1.;
//   return;
// }
// template <> inline void phiBase<2>(double *array, const double x) noexcept {
//   array[0] = (1. - x);
//   array[1] = x;
//   return;
// }
// template <> inline void phiBase<3>(double *array, const double x) noexcept {

//   array[0] = (1. - x) * (1. - x);
//   array[1] = (2. * x) * (1. - x);
//   array[2] = x * x;
//   return;
// }
// template <> inline void phiBase<4>(double *array, const double x) noexcept {

//   array[0] = (1. - x) * (1. - x) * (1. - x);
//   array[1] = 3.0 * x * (1. - x) * (1. - x);
//   array[2] = 3. * x * x * (1. - x);
//   array[3] = x * x * x;
//   return;
// }
// template <> inline void phiBase<5>(double *array, const double x) noexcept {

//   array[0] = (1. - x) * (1. - x) * (1. - x) * (1. - x);
//   array[1] = 4. * x * (1. - x) * (1. - x) * (1. - x);
//   array[2] = 6. * x * x * (1. - x) * (1. - x);
//   array[3] = 4. * x * x * x * (1. - x);
//   array[4] = x * x * x * x;
//   return;
// }

// //************************************************************************
// //********* PHIDER
// //************************************************************************

// template <> inline void phiBaseDer<1>(double *array, const double x) noexcept
// {
//   array[0] = 0.;
//   return;
// }
// template <> inline void phiBaseDer<2>(double *array, const double x) noexcept
// {
//   array[0] = -1.;
//   array[1] = 1;
//   return;
// }
// template <> inline void phiBaseDer<3>(double *array, const double x) noexcept
// {

//   array[0] = 2. * x - 2;
//   array[1] = 2. - 4. * x;
//   array[2] = 2. * x;
//   return;
// }
// template <> inline void phiBaseDer<4>(double *array, const double x) noexcept
// {

//   array[0] = -3 * (1. - x) * (1. - x);
//   array[1] = 3. - 12. * x + 9. * x * x;
//   array[2] = (6. - 9. * x) * x;
//   array[3] = 3. * x * x;
//   return;
// }
// template <> inline void phiBaseDer<5>(double *array, const double x) noexcept
// {

//   array[0] = -4 * (1. - x) * (1. - x) * (1. - x);
//   array[1] = -16. * (1. - x) * (1. - x) * (-0.25 + x);
//   array[2] = x * (12. - 36. * x + 24. * x * x);
//   array[3] = (12 - 16 * x) * x * x;
//   array[4] = 4. * x * x * x;
//   return;
// }

template <int I>
void phi(double *c, double w, double x) noexcept {
  double base[I];
  Spl::evalBrnstnBasis<double, I - 1>(base, x);
  for (int i = 0; i < I; i++) c[i] += w * base[i];
  return;
}

template <int I>
void phi_dx(double *c, double w, double x) noexcept {
  double base[I];
  Spl::evalBrnstnDerBasis<double, I - 1>(base, x);
  for (int i = 0; i < I; i++) c[i] += w * base[i];
  return;
}

template <int I>
void phiphi(double *c, vector2 a) noexcept {
  double X[I], Y[I];
  Spl::evalBrnstnBasis<double, I - 1>(X, a.x);
  Spl::evalBrnstnBasis<double, I - 1>(Y, a.y);

  for (int iy = 0; iy < I; iy++)
    for (int ix = 0; ix < I; ix++) c[iy * I + ix] = X[ix] * Y[iy];

  return;
}

template <int I>
void phiphi_dx(double *c, vector2 a) noexcept {
  double dX[I], Y[I];
  Spl::evalBrnstnDerBasis<double, I - 1>(dX, a.x);
  Spl::evalBrnstnBasis<double, I - 1>(Y, a.y);

  for (int iy = 0; iy < I; iy++)
    for (int ix = 0; ix < I; ix++) c[iy * I + ix] = dX[ix] * Y[iy];

  return;
}

template <int I>
void phiphi_dy(double *c, vector2 a) noexcept {
  double X[I], dY[I];
  Spl::evalBrnstnBasis<double, I - 1>(X, a.x);
  Spl::evalBrnstnDerBasis<double, I - 1>(dY, a.y);

  for (int iy = 0; iy < I; iy++)
    for (int ix = 0; ix < I; ix++) c[iy * I + ix] = X[ix] * dY[iy];
  return;
}

template <int I>
void Phi_times_Phi(double *c, double w, vector2 xi, vector2 eta) noexcept {
  double a[I * I], b[I * I], X[I], Y[I];

  Spl::evalBrnstnBasis<double, I - 1>(X, xi.x);
  Spl::evalBrnstnBasis<double, I - 1>(Y, xi.y);

  for (int iy = 0; iy < I; iy++)
    for (int ix = 0; ix < I; ix++) a[iy * I + ix] = w * X[ix] * Y[iy];

  Spl::evalBrnstnBasis<double, I - 1>(X, eta.x);
  Spl::evalBrnstnBasis<double, I - 1>(Y, eta.y);

  for (int iy = 0; iy < I; iy++)
    for (int ix = 0; ix < I; ix++) b[iy * I + ix] = X[ix] * Y[iy];

  for (int i = 0; i < (I * I); i++)
    for (int j = 0; j < (I * I); j++) c[i * (I * I) + j] += a[i] * b[j];

  return;
}

template <int I>
void VPhi_scal_VPhi(double *c, double weight, vector2 xi, vector2 eta,
                    vector3 dChi_dx_s, vector3 dChi_dy_s, vector3 dChi_dx_t,
                    vector3 dChi_dy_t)
// double *c,
// weight;
// vector2 xi, eta;
// vector3 dChi_dx_s, dChi_dy_s, dChi_dx_t, dChi_dy_t;
{
  constexpr int a_bs = I * I;
  constexpr int a_bs2 = a_bs * a_bs;

  double d[a_bs2];

  memset(d, 0, a_bs2 * sizeof(double));
  Phi_times_Phi<I>(d, weight, xi, eta);

  mydaxpy(a_bs2, vector3_skalp(dChi_dx_s, dChi_dx_t), d, c);
  mydaxpy(a_bs2, vector3_skalp(dChi_dx_s, dChi_dy_t), d, c + a_bs2);
  mydaxpy(a_bs2, vector3_skalp(dChi_dy_s, dChi_dx_t), d, c + 2 * a_bs2);
  mydaxpy(a_bs2, vector3_skalp(dChi_dy_s, dChi_dy_t), d, c + 3 * a_bs2);

  return;
}

/**
 *  \brief         Updates \f$c_{ij}\f$ by \f$\mathrm{weight}\cdot
 *                 \div_{\Gamma}\phi _i(\mathrm{xi})\cdot\div_{\Gamma}\phi
 * _j(\mathrm{eta})\f$.
 *
 *  \author        JD&FW.
 */
template <int I>
void Div_Phi_times_Div_Phi(double *c, double weight, vector2 xi, vector2 eta)
// double *c;
// double weight;
// vector2 xi;
// vector2 eta;
{
  double a_dx[I * I];
  double a_dy[I * I];
  double b_dx[I * I];
  double b_dy[I * I];

  phiphi_dx<I>(a_dx, xi);
  mydscal(I * I, weight, a_dx);
  phiphi_dy<I>(a_dy, xi);
  mydscal(I * I, weight, a_dy);
  phiphi_dx<I>(b_dx, eta);
  phiphi_dy<I>(b_dy, eta);

  for (int i = 0; i < I * I; ++i)
    for (int j = 0; j < I * I; ++j) c[I * I * i + j + 0] += a_dx[i] * b_dx[j];
  for (int i = 0; i < I * I; ++i)
    for (int j = 0; j < I * I; ++j)
      c[I * I * i + j + I * I * I * I] += a_dx[i] * b_dy[j];
  for (int i = 0; i < I * I; ++i)
    for (int j = 0; j < I * I; ++j)
      c[I * I * i + j + 2 * I * I * I * I] += a_dy[i] * b_dx[j];
  for (int i = 0; i < I * I; ++i)
    for (int j = 0; j < I * I; ++j)
      c[I * I * i + j + 3 * I * I * I * I] += a_dy[i] * b_dy[j];

  return;
}

typedef void (*_tmp_phi_type_)(double *, double, double);
typedef void (*_tmp_phi_dx_type_)(double *, double, double);
typedef void (*_tmp_phiphi_type_)(double *, vector2);
typedef void (*_tmp_phi_times_phi_type_)(double *, double, vector2, vector2);
typedef void (*_VPhi_scal_VPhi_type_)(double *, double, vector2, vector2,
                                      vector3, vector3, vector3, vector3);
typedef void (*_Div_Phi_times_Div_Phi_type_)(double *, double, vector2,
                                             vector2);

inline _tmp_phi_type_ select_phi(int deg) {
  switch (deg) {
    case (1):
      return &phi<1>;
    case (2):
      return &phi<2>;
    case (3):
      return &phi<3>;
    case (4):
      return &phi<4>;
    case (5):
      return &phi<5>;
    case (6):
      return &phi<6>;
    case (7):
      return &phi<7>;
    case (8):
      return &phi<8>;
    case (9):
      return &phi<9>;
    case (10):
      return &phi<10>;
    case (11):
      return &phi<11>;
    case (12):
      return &phi<12>;
    case (13):
      return &phi<13>;
    case (14):
      return &phi<14>;
    case (15):
      return &phi<15>;
    case (16):
      return &phi<16>;
    case (17):
      return &phi<17>;
    case (18):
      return &phi<18>;
    default:
      return _tmp_phi_type_(NULL);
  }
}
inline _tmp_phi_dx_type_ select_phi_dx(int deg) {
  switch (deg) {
    case (1):
      return &phi_dx<1>;
    case (2):
      return &phi_dx<2>;
    case (3):
      return &phi_dx<3>;
    case (4):
      return &phi_dx<4>;
    case (5):
      return &phi_dx<5>;
    case (6):
      return &phi_dx<6>;
    case (7):
      return &phi_dx<7>;
    case (8):
      return &phi_dx<8>;
    case (9):
      return &phi_dx<9>;
    case (10):
      return &phi_dx<10>;
    case (11):
      return &phi_dx<11>;
    case (12):
      return &phi_dx<12>;
    case (13):
      return &phi_dx<13>;
    case (14):
      return &phi_dx<14>;
    case (15):
      return &phi_dx<15>;
    case (16):
      return &phi_dx<16>;
    case (17):
      return &phi_dx<17>;
    case (18):
      return &phi_dx<18>;
    default:
      return _tmp_phi_dx_type_(NULL);
  }
}
inline _tmp_phiphi_type_ select_phiphi(int deg) {
  switch (deg) {
    case (1):
      return &phiphi<1>;
    case (2):
      return &phiphi<2>;
    case (3):
      return &phiphi<3>;
    case (4):
      return &phiphi<4>;
    case (5):
      return &phiphi<5>;
    case (6):
      return &phiphi<6>;
    case (7):
      return &phiphi<7>;
    case (8):
      return &phiphi<8>;
    case (9):
      return &phiphi<9>;
    case (10):
      return &phiphi<10>;
    case (11):
      return &phiphi<11>;
    case (12):
      return &phiphi<12>;
    case (13):
      return &phiphi<13>;
    case (14):
      return &phiphi<14>;
    case (15):
      return &phiphi<15>;
    case (16):
      return &phiphi<16>;
    case (17):
      return &phiphi<17>;
    case (18):
      return &phiphi<18>;
    default:
      return _tmp_phiphi_type_(NULL);
  }
}
inline _tmp_phiphi_type_ select_phiphi_dx(int deg) {
  switch (deg) {
    case (1):
      return &phiphi_dx<1>;
    case (2):
      return &phiphi_dx<2>;
    case (3):
      return &phiphi_dx<3>;
    case (4):
      return &phiphi_dx<4>;
    case (5):
      return &phiphi_dx<5>;
    case (6):
      return &phiphi_dx<6>;
    case (7):
      return &phiphi_dx<7>;
    case (8):
      return &phiphi_dx<8>;
    case (9):
      return &phiphi_dx<9>;
    case (10):
      return &phiphi_dx<10>;
    case (11):
      return &phiphi_dx<11>;
    case (12):
      return &phiphi_dx<12>;
    case (13):
      return &phiphi_dx<13>;
    case (14):
      return &phiphi_dx<14>;
    case (15):
      return &phiphi_dx<15>;
    case (16):
      return &phiphi_dx<16>;
    case (17):
      return &phiphi_dx<17>;
    case (18):
      return &phiphi_dx<18>;
    default:
      return _tmp_phiphi_type_(NULL);
  }
}
inline _tmp_phiphi_type_ select_phiphi_dy(int deg) {
  switch (deg) {
    case (1):
      return &phiphi_dy<1>;
    case (2):
      return &phiphi_dy<2>;
    case (3):
      return &phiphi_dy<3>;
    case (4):
      return &phiphi_dy<4>;
    case (5):
      return &phiphi_dy<5>;
    case (6):
      return &phiphi_dy<6>;
    case (7):
      return &phiphi_dy<7>;
    case (8):
      return &phiphi_dy<8>;
    case (9):
      return &phiphi_dy<9>;
    case (10):
      return &phiphi_dy<10>;
    case (11):
      return &phiphi_dy<11>;
    case (12):
      return &phiphi_dy<12>;
    case (13):
      return &phiphi_dy<13>;
    case (14):
      return &phiphi_dy<14>;
    case (15):
      return &phiphi_dy<15>;
    case (16):
      return &phiphi_dy<16>;
    case (17):
      return &phiphi_dy<17>;
    case (18):
      return &phiphi_dy<18>;
    default:
      return _tmp_phiphi_type_(NULL);
  }
}
inline _tmp_phi_times_phi_type_ select_Phi_times_Phi(int deg) {
  switch (deg) {
    case (1):
      return &Phi_times_Phi<1>;
    case (2):
      return &Phi_times_Phi<2>;
    case (3):
      return &Phi_times_Phi<3>;
    case (4):
      return &Phi_times_Phi<4>;
    case (5):
      return &Phi_times_Phi<5>;
    case (6):
      return &Phi_times_Phi<6>;
    case (7):
      return &Phi_times_Phi<7>;
    case (8):
      return &Phi_times_Phi<8>;
    case (9):
      return &Phi_times_Phi<9>;
    case (10):
      return &Phi_times_Phi<10>;
    case (11):
      return &Phi_times_Phi<11>;
    case (12):
      return &Phi_times_Phi<12>;
    case (13):
      return &Phi_times_Phi<13>;
    case (14):
      return &Phi_times_Phi<14>;
    case (15):
      return &Phi_times_Phi<15>;
    case (16):
      return &Phi_times_Phi<16>;
    case (17):
      return &Phi_times_Phi<17>;
    case (18):
      return &Phi_times_Phi<18>;
    default:
      return _tmp_phi_times_phi_type_(NULL);
  }
}

inline _VPhi_scal_VPhi_type_ select_VPhi_scal_bla(int deg) {
  switch (deg) {
    case (1):
      return &VPhi_scal_VPhi<1>;
    case (2):
      return &VPhi_scal_VPhi<2>;
    case (3):
      return &VPhi_scal_VPhi<3>;
    case (4):
      return &VPhi_scal_VPhi<4>;
    case (5):
      return &VPhi_scal_VPhi<5>;
    case (6):
      return &VPhi_scal_VPhi<6>;
    case (7):
      return &VPhi_scal_VPhi<7>;
    case (8):
      return &VPhi_scal_VPhi<8>;
    case (9):
      return &VPhi_scal_VPhi<9>;
    case (10):
      return &VPhi_scal_VPhi<10>;
    case (11):
      return &VPhi_scal_VPhi<11>;
    case (12):
      return &VPhi_scal_VPhi<12>;
    case (13):
      return &VPhi_scal_VPhi<13>;
    case (14):
      return &VPhi_scal_VPhi<14>;
    case (15):
      return &VPhi_scal_VPhi<15>;
    case (16):
      return &VPhi_scal_VPhi<16>;
    case (17):
      return &VPhi_scal_VPhi<17>;
    case (18):
      return &VPhi_scal_VPhi<18>;
    default:
      return _VPhi_scal_VPhi_type_(NULL);
  }
}
inline _Div_Phi_times_Div_Phi_type_ select_Div_Phi_times_Div_Phi(int deg) {
  switch (deg) {
    case (1):
      return &Div_Phi_times_Div_Phi<1>;
    case (2):
      return &Div_Phi_times_Div_Phi<2>;
    case (3):
      return &Div_Phi_times_Div_Phi<3>;
    case (4):
      return &Div_Phi_times_Div_Phi<4>;
    case (5):
      return &Div_Phi_times_Div_Phi<5>;
    case (6):
      return &Div_Phi_times_Div_Phi<6>;
    case (7):
      return &Div_Phi_times_Div_Phi<7>;
    case (8):
      return &Div_Phi_times_Div_Phi<8>;
    case (9):
      return &Div_Phi_times_Div_Phi<9>;
    case (10):
      return &Div_Phi_times_Div_Phi<10>;
    case (11):
      return &Div_Phi_times_Div_Phi<11>;
    case (12):
      return &Div_Phi_times_Div_Phi<12>;
    case (13):
      return &Div_Phi_times_Div_Phi<13>;
    case (14):
      return &Div_Phi_times_Div_Phi<14>;
    case (15):
      return &Div_Phi_times_Div_Phi<15>;
    case (16):
      return &Div_Phi_times_Div_Phi<16>;
    case (17):
      return &Div_Phi_times_Div_Phi<17>;
    case (18):
      return &Div_Phi_times_Div_Phi<18>;
    default:
      return _Div_Phi_times_Div_Phi_type_(NULL);
  }
}
}  // namespace Bembel
#endif
