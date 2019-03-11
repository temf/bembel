// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "spline/bernstein.h"

namespace Spl {

/**
 *	@brief Routines for basis function evaluation
 *
 *	Here we setup functions for the evaluation of the basis. PEQ -> pointer
 *equals, meaning += instead of =  is used in value assignment. Different
 *formats are supported, see type signatures
 */
std::vector<double> evalBrnstn(int deg, double *coef,
                               const std::vector<double> &x);

double evalBrnstn(int deg, double *coef, double x);

void evalBrnstnBasis(int deg, double *out, double x);

std::vector<double> evalBrnstnBasis(int deg, double x);

void evalBrnstnBasisPEQ(int deg, double *out, double x);

std::vector<double> evalBrnstnDer(int deg, double *coef,
                                  const std::vector<double> &x);

double evalBrnstnDer(int deg, double *coef, double x);

void evalBrnstnDerBasis(int deg, double *out, double x);

void evalBrnstnDerBasisPEQ(int deg, double *out, double x);

std::vector<double> evalBrnstnDerBasis(int deg, double x);
}  // namespace Spl