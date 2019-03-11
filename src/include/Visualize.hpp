// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_VISUALIZE_
#define __BEMBEL_VISUALIZE_

#include "Conversions.hpp"
#include "Discretization.hpp"
#include "Eigen/Dense"
#include "Spline.hpp"

namespace Bembel {
namespace Vis {

void plot_geometry_only(const std::vector<Spl::Patch> &geoms, int num,
                        const char *);
void plot_rho_laplace(Bembel::Discretization<Bembel::LaplaceSingle> &myDisc,
                      const Eigen::VectorXd &rho, int num, const char *name);

void plot_rho_helmholtz(Bembel::Discretization<Bembel::HelmholtzSingle> &myDisc,
                        const Eigen::VectorXcd &rho, int num, const char *name);

void plot_rho_maxwell(Bembel::Discretization<Bembel::MaxwellSingle> &myDisc,
                      const Eigen::VectorXcd &rho, int num, const char *name);

template <typename T>
inline void plotDiscretizationToVTK(Bembel::Discretization<T> &myDisc,
                                    std::string name = "geometry.vtk",
                                    int num = 5) {
  plot_geometry_only(myDisc.get_plain_patchdata(), num, name.c_str());
  return;
}

inline void plotDiscretizationToVTK(
    Bembel::Discretization<Bembel::LaplaceSingle> &myDisc,
    const Eigen::VectorXd &rho, std::string name = "laplace.vtk", int num = 5) {
  plot_rho_laplace(myDisc, rho, num, name.c_str());

  return;
}

inline void plotDiscretizationToVTK(
    Bembel::Discretization<Bembel::HelmholtzSingle> &myDisc,
    const Eigen::VectorXcd &rho, std::string name = "helmholtz.vtk",
    int num = 5) {
  plot_rho_helmholtz(myDisc, rho, num, name.c_str());

  return;
}

inline void plotDiscretizationToVTK(
    Bembel::Discretization<Bembel::MaxwellSingle> &myDisc,
    const Eigen::VectorXcd &rho, std::string name = "maxwell.vtk",
    int num = 5) {
  plot_rho_maxwell(myDisc, rho, num, name.c_str());
  return;
}

}  // namespace Vis
}  // namespace Bembel

#endif
