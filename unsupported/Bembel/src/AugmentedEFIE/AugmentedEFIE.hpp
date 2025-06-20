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
#ifndef BEMBEL_SRC_AUGMENTEDEFIE_AUGMENTEDEFIE_HPP_
#define BEMBEL_SRC_AUGMENTEDEFIE_AUGMENTEDEFIE_HPP_

namespace Bembel {
/**
 * \ingroup AugmentedEFIE
 * \brief This class handles the assembly and storage of the system matrix.
 */
template <typename MatrixFormat, typename Derived>
class AugmentedEFIE {
 public:
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  AugmentedEFIE() {}
  AugmentedEFIE(const AnsatzSpace<Derived> &ansatz_space_vector,
                const Geometry &geometry) {
    init_AugmentedEFIE(ansatz_space_vector, geometry);
  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_AugmentedEFIE
  //////////////////////////////////////////////////////////////////////////////
  void init_AugmentedEFIE(const AnsatzSpace<Derived> &ansatz_space_vector,
                          const Geometry &geometry) {
    ansatz_space_vector_ = ansatz_space_vector;

    // Build ansatz spaces
    ansatz_space_mass_ = AnsatzSpace<MassMatrixScalarDisc>(
        geometry, ansatz_space_vector.get_refinement_level(),
        ansatz_space_vector.get_polynomial_degree() - 1);

    ansatz_space_scalar_ = AnsatzSpace<HelmholtzSingleLayerOperator>(
        geometry, ansatz_space_vector.get_refinement_level(),
        ansatz_space_vector.get_polynomial_degree() - 1);

    dofs_scalar_ = ansatz_space_scalar_.get_number_of_dofs();
    dofs_vector_ = ansatz_space_vector_.get_number_of_dofs();
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    compute
  //////////////////////////////////////////////////////////////////////////////
  /**
   *  \brief Assembles the system matrix without voltage source excitation.
   */
  inline void compute() {
    system_matrix_ = Eigen::MatrixXcd(dofs_scalar_ + dofs_vector_,
                                      dofs_scalar_ + dofs_vector_);

    AugmentedEFIEAssembler<MatrixFormat>::compute(
        &system_matrix_, ansatz_space_mass_, ansatz_space_scalar_,
        ansatz_space_vector_, wavenumber_);
    return;
  }
  /**
   *  \brief Assembles the system matrix with voltage source excitation.
   */
  inline void compute(VoltageSource source) {
    const int elements_per_edge =
        (1 << ansatz_space_vector_.get_refinement_level());
    // This number contains all degrees of freedom on the edge for the
    // excitation and three extra ones dedicated to the potential of the voltage
    // source.
    const int dofs_excitation =
        3 +
        (source.positive_.size() + source.negative_.size()) * elements_per_edge;
    system_matrix_ =
        Eigen::MatrixXcd::Zero(dofs_scalar_ + dofs_vector_ + dofs_excitation,
                               dofs_scalar_ + dofs_vector_ + dofs_excitation);
    excitation_ = Eigen::VectorXcd::Zero(dofs_scalar_);

    AugmentedEFIEAssembler<MatrixFormat>::compute_matrix_and_excitation(
        &system_matrix_, &excitation_, ansatz_space_mass_, ansatz_space_scalar_,
        ansatz_space_vector_, wavenumber_, source);
    return;
  }
  void stabilize() {
    std::function<std::complex<double>(Eigen::Vector3d)> one =
        [](Eigen::Vector3d in) { return std::complex<double>(1., 0.); };
    DiscreteLinearForm<DirichletTrace<std::complex<double>>,
                       HelmholtzSingleLayerOperator>
        const_lf(ansatz_space_scalar_);
    const_lf.get_linear_form().set_function(one);
    const_lf.compute();

    system_matrix_.bottomRows(dofs_scalar_) *= omega_ *
                                               std::complex<double>(0., 1.) *
                                               Constants::mu0 * Constants::eps0;
    system_matrix_.leftCols(dofs_vector_) *=
        1. / (omega_ * std::complex<double>(0., 1.) * Constants::mu0);

    std::complex<double> eigenvalues_average =
        system_matrix_.trace() / ((double)system_matrix_.cols());

    Eigen::VectorXcd a = Eigen::VectorXcd::Zero(system_matrix_.cols());
    a.bottomRows(dofs_scalar_) = const_lf.get_discrete_linear_form();

    system_matrix_ = system_matrix_ - eigenvalues_average * a * a.transpose();
  }
  //////////////////////////////////////////////////////////////////////////////
  //    setters
  //////////////////////////////////////////////////////////////////////////////
  /**
   *  \brief Set wave number
   */
  void set_wavenumber(std::complex<double> wavenumber) {
    wavenumber_ = wavenumber;
  }
  void set_omega(double omega) { omega_ = omega; }
  //////////////////////////////////////////////////////////////////////////////
  //    getters
  //////////////////////////////////////////////////////////////////////////////
  const MatrixFormat &get_system_matrix() const { return system_matrix_; }
  /**
   *  \brief Return wave number
   */
  const std::complex<double> get_wavenumber() { return wavenumber_; }
  /**
   *  \brief Return number of degrees of freedom for the current.
   */
  const int get_dofs_vector() { return dofs_vector_; }
  /**
   *  \brief Return number of degrees of freedom for the potential
   */
  const int get_dofs_scalar() { return dofs_scalar_; }
  /**
   *  \brief Debug output for the excitation.
   */
  const Eigen::VectorXcd get_excitation() { return excitation_; }
  /**
   *  \brief Debug output for the excitation.
   */
  const AnsatzSpace<HelmholtzSingleLayerOperator> get_ansatz_space_scalar() {
    return ansatz_space_scalar_;
  }
  const AnsatzSpace<MassMatrixScalarDisc> get_ansatz_space_mass() {
    return ansatz_space_mass_;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    private member variables
  //////////////////////////////////////////////////////////////////////////////
 private:
  MatrixFormat system_matrix_;
  Eigen::VectorXcd excitation_;

  AnsatzSpace<MassMatrixScalarDisc> ansatz_space_mass_;
  AnsatzSpace<HelmholtzSingleLayerOperator> ansatz_space_scalar_;
  AnsatzSpace<InductanceMatrix> ansatz_space_vector_;

  int dofs_scalar_;
  int dofs_vector_;
  double omega_;
  std::complex<double> wavenumber_;
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_AUGMENTEDEFIE_AUGMENTEDEFIE_HPP_
