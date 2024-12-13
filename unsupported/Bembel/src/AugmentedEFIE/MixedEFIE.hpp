// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2024 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_SRC_AUGMENTEDEFIE_MIXEDEFIE_HPP_
#define BEMBEL_SRC_AUGMENTEDEFIE_MIXEDEFIE_HPP_

namespace Bembel {
// forward declaration of class MixedEFIE in order to define
// traits

/**
 * \ingroup AugmentedEFIE
 */
template <typename LinOpVector, typename LinOpScalar>
class MixedEFIE {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  MixedEFIE() {}
  MixedEFIE(AnsatzSpace<LinOpVector> ansatz_space_vector,
            AnsatzSpace<LinOpScalar> ansatz_space_scalar) {
    init_MixedEFIE(ansatz_space_vector, ansatz_space_scalar);
  }

  void init_MixedEFIE(AnsatzSpace<LinOpVector> ansatz_space_vector,
                      AnsatzSpace<LinOpScalar> ansatz_space_scalar) {
    ansatz_space_vector_ = ansatz_space_vector;
    ansatz_space_scalar_ = ansatz_space_scalar;
  }

  Eigen::Matrix<std::complex<double>, -1, 3> evaluate(
      const Eigen::Matrix<double, Eigen::Dynamic, 3> &points) {
    // Define the two parts of the EFIE each as Potential -> new classes
    DiscretePotential<VectorPotential<LinOpVector>, LinOpVector> A(
        ansatz_space_vector_);
    DiscretePotential<GradScalarPotential<LinOpScalar>, LinOpScalar> grad_phi(
        ansatz_space_scalar_);

    // Set wave number
    A.get_potential().set_wavenumber(wavenumber_);
    grad_phi.get_potential().set_wavenumber(wavenumber_);

    // Set cauchy data
    A.set_cauchy_data(current_);
    grad_phi.set_cauchy_data(charges_);

    // Compute each contribution to the potential
    Eigen::Matrix<std::complex<double>, -1, 3> potential(points.rows(), 3);
    potential.setZero();
    potential = -std::complex<double>(0., 1.) * omega_ * Constants::mu0 *
                A.evaluate(points);
    potential -= grad_phi.evaluate(points) / Constants::eps0;

    return potential;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    setters
  //////////////////////////////////////////////////////////////////////////////
  void set_omega(double omega) { omega_ = omega; }
  void set_wavenumber(std::complex<double> wavenumber) {
    wavenumber_ = wavenumber;
  }
  void set_current(Eigen::VectorXcd current) { current_ = current; }
  void set_charges(Eigen::VectorXcd charges) { charges_ = charges; }
  //////////////////////////////////////////////////////////////////////////////
  //    getters
  //////////////////////////////////////////////////////////////////////////////
  double get_omega() { return omega_; }
  std::complex<double> get_wavenumber() { return wavenumber_; }
  Eigen::VectorXcd get_current() { return current_; }
  Eigen::VectorXcd get_charges() { return charges_; }

 private:
  AnsatzSpace<LinOpVector> ansatz_space_vector_;
  AnsatzSpace<LinOpScalar> ansatz_space_scalar_;
  Eigen::VectorXcd current_;
  Eigen::VectorXcd charges_;
  double omega_;
  std::complex<double> wavenumber_;
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_AUGMENTEDEFIE_MIXEDEFIE_HPP_
