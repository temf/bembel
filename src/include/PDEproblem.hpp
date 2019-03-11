// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_PDEPROBLEM__
#define __BEMBEL_PDEPROBLEM__

#include "pdeproblem.h"

namespace Bembel {
/*
 * Here are the classes which contain the basics about the PDE.
 *
 */

class LaplaceSingle;
class HelmholtzSingle;
class MaxwellSingle;

template <typename Derived>
struct PDEproblemTraits {};

// specialize typedefs for different PDEs
template <>
struct PDEproblemTraits<LaplaceSingle> {
  typedef Eigen::VectorXd EigenType;
  typedef Eigen::VectorXd::Scalar Scalar;
};

template <>
struct PDEproblemTraits<HelmholtzSingle> {
  typedef Eigen::VectorXcd EigenType;
  typedef Eigen::VectorXcd::Scalar Scalar;
};

template <>
struct PDEproblemTraits<MaxwellSingle> {
  typedef Eigen::VectorXcd EigenType;
  typedef Eigen::VectorXcd::Scalar Scalar;
};

template <typename Derived>
struct PDEproblemBase {
  // Constructors
  PDEproblemBase(){};

  typename PDEproblemTraits<Derived>::EigenType MatVec(
      void *H, const typename PDEproblemTraits<Derived>::EigenType &x) const {
    return static_cast<const Derived *>(this)->MatVecImplementation(H, x);
  };
  pdeproblem _pde;
};

class LaplaceSingle : public PDEproblemBase<LaplaceSingle> {
 public:
  LaplaceSingle(const PDEproblemBase<LaplaceSingle> &other) {
    _pde = other._pde;
  }
  LaplaceSingle(PDEproblemBase<LaplaceSingle> &&other) {
    _pde = std::move(other._pde);
  }
  LaplaceSingle &operator=(PDEproblemBase<LaplaceSingle> other) {
    std::swap(_pde, other._pde);
    return *this;
  }
  Eigen::VectorXd MatVecImplementation(void *, const Eigen::VectorXd &) const;
  LaplaceSingle(void);
  void update_wavenumber(std::complex<double> in) {
    assert(false && "Wavenumber not applicable to LaplaceSingle");
  }
  void update_eta_by_wavenumber(double &in) { return; }
};

class HelmholtzSingle : public PDEproblemBase<HelmholtzSingle> {
 public:
  HelmholtzSingle(const PDEproblemBase<HelmholtzSingle> &other) {
    _pde = other._pde;
  }
  HelmholtzSingle(PDEproblemBase<HelmholtzSingle> &&other) {
    _pde = std::move(other._pde);
  }
  HelmholtzSingle &operator=(PDEproblemBase<HelmholtzSingle> other) {
    std::swap(_pde, other._pde);
    return *this;
  }
  Eigen::VectorXcd MatVecImplementation(void *, const Eigen::VectorXcd &) const;
  HelmholtzSingle(std::complex<double> in = std::complex<double>(1, 0));
  void update_wavenumber(std::complex<double> in) {
    _pde.kappa[0] = in.real();
    _pde.kappa[1] = in.imag();
  }

  std::complex<double> get_wavenumber() {
    return std::complex<double>(_pde.kappa[0], _pde.kappa[1]);
  }
  void update_eta_by_wavenumber(double &in) {
    in /= std::abs(get_wavenumber());
  }
};

class MaxwellSingle : public PDEproblemBase<MaxwellSingle> {
 public:
  MaxwellSingle(const PDEproblemBase<MaxwellSingle> &other) {
    _pde = other._pde;
  }
  MaxwellSingle(PDEproblemBase<MaxwellSingle> &&other) {
    _pde = std::move(other._pde);
  }
  MaxwellSingle &operator=(PDEproblemBase<MaxwellSingle> other) {
    std::swap(_pde, other._pde);
    return *this;
  }
  Eigen::VectorXcd MatVecImplementation(void *, const Eigen::VectorXcd &) const;
  MaxwellSingle(std::complex<double> in = std::complex<double>(1, 0));
  void update_wavenumber(std::complex<double> in) {
    _pde.kappa[0] = in.real();
    _pde.kappa[1] = in.imag();
  }
  std::complex<double> get_wavenumber() {
    return std::complex<double>(_pde.kappa[0], _pde.kappa[1]);
  }
  void update_eta_by_wavenumber(double &in) {
    in /= std::abs(get_wavenumber());
  }
};
}  // namespace Bembel
#endif
