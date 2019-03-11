// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_DISCRETIZATION__
#define __BEMBEL_DISCRETIZATION__

#include "Geometry.hpp"
#include "Mesh.hpp"
#include "PDEproblem.hpp"
#include "discretization.h"

namespace Bembel {
template <typename Derived>
class Discretization {
 public:
  Discretization() {}
  inline Discretization(const Geometry& geo, const PDEproblemBase<Derived>& pde,
                        int deg, int kntrep, int M) {
    _geo = geo;
    _msh.init_mesh(_geo, M);
    _pde = pde;
    _disc =
        get_discretization_NB(deg, kntrep, &(_pde._pde), &(_msh.get_mesh()));
  }
  ~Discretization() { free_discretization(&_disc); }
  void init_Discretization(const Geometry& geo,
                           const PDEproblemBase<Derived>& pde, int deg,
                           int kntrep, int M) {
    _geo = geo;
    _msh.init_mesh(_geo, M);
    _pde = pde;
    _disc =
        get_discretization_NB(deg, kntrep, &(_pde._pde), &(_msh.get_mesh()));
  }
  discretization& get_disc() { return _disc; }
  Derived& get_pde() { return _pde; }
  std::vector<Spl::Patch>& get_plain_patchdata() { return _geo.get_geometry(); }
  void update_wavenumber(std::complex<double> in) {
    _pde.update_wavenumber(in);
    return;
  }
  int get_num_dofs() { return _disc.real_na; }

 private:
  discretization _disc;
  Geometry _geo;
  Mesh _msh;
  Derived _pde;
};

}  // namespace Bembel
#endif
