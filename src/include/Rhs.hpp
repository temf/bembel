// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_RHS_
#define __BEMBEL_RHS_
#include <Eigen/Dense>
#include <functional>
#include "BEMRHS.h"
#include "Conversions.hpp"
#include "Discretization.hpp"
#include "PDEproblem.hpp"
#include "Spline.hpp"
#include "myvector.h"

namespace Bembel {
namespace Rhs {
Eigen::VectorXd computeRhs(Discretization<LaplaceSingle> &ddisc,
                           std::function<double(Eigen::Vector3d)> fun);

Eigen::VectorXcd computeRhs(
    Discretization<HelmholtzSingle> &ddisc,
    std::function<std::complex<double>(Eigen::Vector3d, std::complex<double>)>
        fun);

Eigen::VectorXcd computeRhs(
    Discretization<MaxwellSingle> &ddisc,
    std::function<Eigen::Vector3cd(Eigen::Vector3d, std::complex<double>)> fun);
}  // namespace Rhs
}  // namespace Bembel
#endif
