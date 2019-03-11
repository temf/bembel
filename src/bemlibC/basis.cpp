// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "spline/bernstein.h"

namespace Spl {

/////// Nasty stuff
typedef double (*fun_bernstein_type_one)(double *, double);
typedef std::vector<double> (*fun_bernstein_type_two)(
    double *, const std::vector<double> &);
typedef void (*fun_bernstein_type_three)(double *, double);

constexpr std::array<fun_bernstein_type_one, 19>
fun_bernstein_type_one_global_routines() {
  return {
      &evalBrnstn<double, 0>,  &evalBrnstn<double, 1>,  &evalBrnstn<double, 2>,
      &evalBrnstn<double, 3>,  &evalBrnstn<double, 4>,  &evalBrnstn<double, 5>,
      &evalBrnstn<double, 6>,  &evalBrnstn<double, 7>,  &evalBrnstn<double, 8>,
      &evalBrnstn<double, 9>,  &evalBrnstn<double, 10>, &evalBrnstn<double, 11>,
      &evalBrnstn<double, 12>, &evalBrnstn<double, 13>, &evalBrnstn<double, 14>,
      &evalBrnstn<double, 15>, &evalBrnstn<double, 16>, &evalBrnstn<double, 17>,
      &evalBrnstn<double, 18>};
}

constexpr std::array<fun_bernstein_type_one, 19>
fun_bernstein_type_one_global_routines_der() {
  return {&evalBrnstnDer<double, 0>,  &evalBrnstnDer<double, 1>,
          &evalBrnstnDer<double, 2>,  &evalBrnstnDer<double, 3>,
          &evalBrnstnDer<double, 4>,  &evalBrnstnDer<double, 5>,
          &evalBrnstnDer<double, 6>,  &evalBrnstnDer<double, 7>,
          &evalBrnstnDer<double, 8>,  &evalBrnstnDer<double, 9>,
          &evalBrnstnDer<double, 10>, &evalBrnstnDer<double, 11>,
          &evalBrnstnDer<double, 12>, &evalBrnstnDer<double, 13>,
          &evalBrnstnDer<double, 14>, &evalBrnstnDer<double, 15>,
          &evalBrnstnDer<double, 16>, &evalBrnstnDer<double, 17>,
          &evalBrnstnDer<double, 18>};
}

constexpr std::array<fun_bernstein_type_two, 19>
fun_bernstein_type_two_global_routines() {
  return {
      &evalBrnstn<double, 0>,  &evalBrnstn<double, 1>,  &evalBrnstn<double, 2>,
      &evalBrnstn<double, 3>,  &evalBrnstn<double, 4>,  &evalBrnstn<double, 5>,
      &evalBrnstn<double, 6>,  &evalBrnstn<double, 7>,  &evalBrnstn<double, 8>,
      &evalBrnstn<double, 9>,  &evalBrnstn<double, 10>, &evalBrnstn<double, 11>,
      &evalBrnstn<double, 12>, &evalBrnstn<double, 13>, &evalBrnstn<double, 14>,
      &evalBrnstn<double, 15>, &evalBrnstn<double, 16>, &evalBrnstn<double, 17>,
      &evalBrnstn<double, 18>};
}

constexpr std::array<fun_bernstein_type_two, 19>
fun_bernstein_type_two_global_routines_der() {
  return {&evalBrnstnDer<double, 0>,  &evalBrnstnDer<double, 1>,
          &evalBrnstnDer<double, 2>,  &evalBrnstnDer<double, 3>,
          &evalBrnstnDer<double, 4>,  &evalBrnstnDer<double, 5>,
          &evalBrnstnDer<double, 6>,  &evalBrnstnDer<double, 7>,
          &evalBrnstnDer<double, 8>,  &evalBrnstnDer<double, 9>,
          &evalBrnstnDer<double, 10>, &evalBrnstnDer<double, 11>,
          &evalBrnstnDer<double, 12>, &evalBrnstnDer<double, 13>,
          &evalBrnstnDer<double, 14>, &evalBrnstnDer<double, 15>,
          &evalBrnstnDer<double, 16>, &evalBrnstnDer<double, 17>,
          &evalBrnstnDer<double, 18>};
}

constexpr std::array<fun_bernstein_type_three, 19>
fun_bernstein_type_three_global_routines() {
  return {&evalBrnstnBasis<double, 0>,  &evalBrnstnBasis<double, 1>,
          &evalBrnstnBasis<double, 2>,  &evalBrnstnBasis<double, 3>,
          &evalBrnstnBasis<double, 4>,  &evalBrnstnBasis<double, 5>,
          &evalBrnstnBasis<double, 6>,  &evalBrnstnBasis<double, 7>,
          &evalBrnstnBasis<double, 8>,  &evalBrnstnBasis<double, 9>,
          &evalBrnstnBasis<double, 10>, &evalBrnstnBasis<double, 11>,
          &evalBrnstnBasis<double, 12>, &evalBrnstnBasis<double, 13>,
          &evalBrnstnBasis<double, 14>, &evalBrnstnBasis<double, 15>,
          &evalBrnstnBasis<double, 16>, &evalBrnstnBasis<double, 17>,
          &evalBrnstnBasis<double, 18>};
}

constexpr std::array<fun_bernstein_type_three, 19>
fun_bernstein_type_three_global_routines_der() {
  return {&evalBrnstnDerBasis<double, 0>,  &evalBrnstnDerBasis<double, 1>,
          &evalBrnstnDerBasis<double, 2>,  &evalBrnstnDerBasis<double, 3>,
          &evalBrnstnDerBasis<double, 4>,  &evalBrnstnDerBasis<double, 5>,
          &evalBrnstnDerBasis<double, 6>,  &evalBrnstnDerBasis<double, 7>,
          &evalBrnstnDerBasis<double, 8>,  &evalBrnstnDerBasis<double, 9>,
          &evalBrnstnDerBasis<double, 10>, &evalBrnstnDerBasis<double, 11>,
          &evalBrnstnDerBasis<double, 12>, &evalBrnstnDerBasis<double, 13>,
          &evalBrnstnDerBasis<double, 14>, &evalBrnstnDerBasis<double, 15>,
          &evalBrnstnDerBasis<double, 16>, &evalBrnstnDerBasis<double, 17>,
          &evalBrnstnDerBasis<double, 18>};
}

constexpr std::array<fun_bernstein_type_three, 19>
fun_bernstein_type_three_global_routines_PEQ() {
  return {&evalBrnstnBasisPEQ<double, 0>,  &evalBrnstnBasisPEQ<double, 1>,
          &evalBrnstnBasisPEQ<double, 2>,  &evalBrnstnBasisPEQ<double, 3>,
          &evalBrnstnBasisPEQ<double, 4>,  &evalBrnstnBasisPEQ<double, 5>,
          &evalBrnstnBasisPEQ<double, 6>,  &evalBrnstnBasisPEQ<double, 7>,
          &evalBrnstnBasisPEQ<double, 8>,  &evalBrnstnBasisPEQ<double, 9>,
          &evalBrnstnBasisPEQ<double, 10>, &evalBrnstnBasisPEQ<double, 11>,
          &evalBrnstnBasisPEQ<double, 12>, &evalBrnstnBasisPEQ<double, 13>,
          &evalBrnstnBasisPEQ<double, 14>, &evalBrnstnBasisPEQ<double, 15>,
          &evalBrnstnBasisPEQ<double, 16>, &evalBrnstnBasisPEQ<double, 17>,
          &evalBrnstnBasisPEQ<double, 18>};
}

constexpr std::array<fun_bernstein_type_three, 19>
fun_bernstein_type_three_global_routines_der_PEQ() {
  return {
      &evalBrnstnDerBasisPEQ<double, 0>,  &evalBrnstnDerBasisPEQ<double, 1>,
      &evalBrnstnDerBasisPEQ<double, 2>,  &evalBrnstnDerBasisPEQ<double, 3>,
      &evalBrnstnDerBasisPEQ<double, 4>,  &evalBrnstnDerBasisPEQ<double, 5>,
      &evalBrnstnDerBasisPEQ<double, 6>,  &evalBrnstnDerBasisPEQ<double, 7>,
      &evalBrnstnDerBasisPEQ<double, 8>,  &evalBrnstnDerBasisPEQ<double, 9>,
      &evalBrnstnDerBasisPEQ<double, 10>, &evalBrnstnDerBasisPEQ<double, 11>,
      &evalBrnstnDerBasisPEQ<double, 12>, &evalBrnstnDerBasisPEQ<double, 13>,
      &evalBrnstnDerBasisPEQ<double, 14>, &evalBrnstnDerBasisPEQ<double, 15>,
      &evalBrnstnDerBasisPEQ<double, 16>, &evalBrnstnDerBasisPEQ<double, 17>,
      &evalBrnstnDerBasisPEQ<double, 18>};
}

std::vector<double> evalBrnstn(int deg, double *coef,
                               const std::vector<double> &x) {
  return fun_bernstein_type_two_global_routines()[deg](coef, x);
}

double evalBrnstn(int deg, double *coef, double x) {
  return fun_bernstein_type_one_global_routines()[deg](coef, x);
}

void evalBrnstnBasis(int deg, double *out, double x) {
  fun_bernstein_type_three_global_routines()[deg](out, x);
  return;
}

std::vector<double> evalBrnstnBasis(int deg, double x) {
  std::vector<double> out(deg + 1);
  evalBrnstnBasis(deg, out.data(), x);
  return out;
}

void evalBrnstnBasisPEQ(int deg, double *out, double x) {
  fun_bernstein_type_three_global_routines_PEQ()[deg](out, x);
  return;
}

std::vector<double> evalBrnstnDer(int deg, double *coef,
                                  const std::vector<double> &x) {
  return fun_bernstein_type_two_global_routines_der()[deg](coef, x);
}

double evalBrnstnDer(int deg, double *coef, double x) {
  return fun_bernstein_type_one_global_routines_der()[deg](coef, x);
}

void evalBrnstnDerBasis(int deg, double *out, double x) {
  fun_bernstein_type_three_global_routines_der()[deg](out, x);
  return;
}

void evalBrnstnDerBasisPEQ(int deg, double *out, double x) {
  fun_bernstein_type_three_global_routines_der_PEQ()[deg](out, x);
  return;
}

std::vector<double> evalBrnstnDerBasis(int deg, double x) {
  std::vector<double> out(deg + 1);
  evalBrnstnDerBasis(deg, out.data(), x);
  return out;
}
}  // namespace Spl