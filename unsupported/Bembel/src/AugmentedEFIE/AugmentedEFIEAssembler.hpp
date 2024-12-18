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
#ifndef BEMBEL_SRC_AUGMENTEDEFIE_AUGMENTEDEFIEASSEMBLER_HPP_
#define BEMBEL_SRC_AUGMENTEDEFIE_AUGMENTEDEFIEASSEMBLER_HPP_

namespace Bembel {
template <typename MatrixFormat>
struct AugmentedEFIEAssembler {};
/**
 * \ingroup AugmentedEFIE
 * \brief This struct handles the routines for assembling the system matrix.
 */
template <>
struct AugmentedEFIEAssembler<Eigen::MatrixXcd> {
 public:
  AugmentedEFIEAssembler() {}
  /**
   * This routine assembles the system matrix for solving the scattering
   * problem.
   */
  static std::unique_ptr<Eigen::MatrixXcd> compute(
      Eigen::MatrixXcd *system_matrix,
      const AnsatzSpace<MassMatrixScalarDisc> &ansatz_space_mass,
      const AnsatzSpace<HelmholtzSingleLayerOperator> &ansatz_space_scalar,
      const AnsatzSpace<InductanceMatrix> &ansatz_space_vector,
      const std::complex<double> wavenumber) {
    const int dofs_scalar = ansatz_space_scalar.get_number_of_dofs();
    const int dofs_vector = ansatz_space_vector.get_number_of_dofs();

    assembleZ(system_matrix, ansatz_space_vector, wavenumber);

    IncidenceMatrix<InductanceMatrix, HelmholtzSingleLayerOperator>
        incidence_matrix(ansatz_space_vector, ansatz_space_scalar);
    incidence_matrix.compute();

    assembleS(system_matrix, incidence_matrix, dofs_vector, dofs_scalar);

    DiscreteLocalOperator<MassMatrixScalarDisc> mass_matrix(ansatz_space_mass);
    mass_matrix.compute();

    assembleM(system_matrix, mass_matrix, wavenumber, dofs_vector, dofs_scalar);

    DiscreteOperator<Eigen::MatrixXcd, HelmholtzSingleLayerOperator>
        capacitive_operator(ansatz_space_scalar);
    capacitive_operator.get_linear_operator().set_wavenumber(wavenumber);
    capacitive_operator.compute();

    std::unique_ptr<Eigen::MatrixXcd> P_times_mass_inverse(
        new Eigen::MatrixXcd);
    *P_times_mass_inverse =
        capacitive_operator.get_discrete_operator() / Constants::eps0 *
        Eigen::MatrixXcd(mass_matrix.get_discrete_operator()).inverse();

    assembleP(system_matrix, ansatz_space_scalar, P_times_mass_inverse,
              incidence_matrix, wavenumber, dofs_vector, dofs_scalar);

    return P_times_mass_inverse;
  }
  /**
   *  This routine assembles the incident matrices necessary for the voltage
   *  source excitation. It is necessary to run the code of the compute routine
   *  above because the inverse of the mass matrix is utilized.
   *
   *  The variable excitation is used for debug output only.
   */
  static void compute_matrix_and_excitation(
      Eigen::MatrixXcd *system_matrix, Eigen::VectorXcd *excitation,
      const AnsatzSpace<MassMatrixScalarDisc> &ansatz_space_mass,
      const AnsatzSpace<HelmholtzSingleLayerOperator> &ansatz_space_scalar,
      const AnsatzSpace<InductanceMatrix> &ansatz_space_vector,
      const std::complex<double> wavenumber, const VoltageSource &source) {
    assert(ansatz_space_scalar.get_polynomial_degree() == 0 &&
           "This works only for lowest order!");
    const int refinement_level = ansatz_space_scalar.get_refinement_level();
    const int dofs_scalar = ansatz_space_scalar.get_number_of_dofs();
    const int dofs_vector = ansatz_space_vector.get_number_of_dofs();

    std::unique_ptr<Eigen::MatrixXcd> P_times_mass_inverse =
        compute(system_matrix, ansatz_space_mass, ansatz_space_scalar,
                ansatz_space_vector, wavenumber);

    assembleV(system_matrix, excitation, ansatz_space_scalar,
              P_times_mass_inverse, source, dofs_vector,
              ansatz_space_scalar.get_number_of_dofs());
    assembleA(system_matrix, source, refinement_level, dofs_vector,
              dofs_scalar);
    return;
  }

 private:
  /**
   *  \brief Assemble the top left block of the system matrix j*omega*mu*L.
   *
   *  Currently there are no finitely resistive materials realized.
   */
  static void assembleZ(
      Eigen::MatrixXcd *system_matrix,
      const AnsatzSpace<InductanceMatrix> &ansatz_space_vector,
      const std::complex<double> wavenumber) {
    DiscreteOperator<Eigen::MatrixXcd, InductanceMatrix> disc_op(
        ansatz_space_vector);
    disc_op.get_linear_operator().set_wavenumber(wavenumber);
    disc_op.compute();

    const int dofs_vector = ansatz_space_vector.get_number_of_dofs();
    std::complex<double> omega =
        wavenumber / std::sqrt(Constants::mu0 * Constants::eps0);
    system_matrix->block(0, 0, dofs_vector, dofs_vector) =
        Constants::mu0 * omega * std::complex<double>(0., 1.) *
        disc_op.get_discrete_operator();
    return;
  }
  /**
   *  \brief Assemble the top right block of the system matrix with the incident
   *  matrix.
   */
  static void assembleS(
      Eigen::MatrixXcd *system_matrix,
      IncidenceMatrix<InductanceMatrix, HelmholtzSingleLayerOperator>
          incidence_matrix,
      const int dofs_vector, const int dofs_scalar) {
    system_matrix->block(0, dofs_vector, dofs_vector, dofs_scalar) =
        incidence_matrix.get_incidence_matrix();
    return;
  }
  /**
   *  \brief Assemble the top right block of the system matrix with the incident
   *  matrix.
   */
  static void assembleM(Eigen::MatrixXcd *system_matrix,
                        DiscreteLocalOperator<MassMatrixScalarDisc> mass_matrix,
                        const std::complex<double> wavenumber,
                        const int dofs_vector, const int dofs_scalar) {
    std::complex<double> omega =
        wavenumber / std::sqrt(Constants::mu0 * Constants::eps0);
    system_matrix->block(dofs_vector, dofs_vector, dofs_scalar, dofs_scalar) =
        -omega * std::complex<double>(0., 1.) *
        mass_matrix.get_discrete_operator();
    return;
  }
  /**
   *  \brief Assemble the bottom left block of the system matrix with
   *  P * M^-1 * S^T.
   */
  static void assembleP(
      Eigen::MatrixXcd *system_matrix,
      const AnsatzSpace<HelmholtzSingleLayerOperator> &ansatz_space_scalar,
      const std::unique_ptr<Eigen::MatrixXcd> &P_times_mass_inverse,
      const IncidenceMatrix<InductanceMatrix, HelmholtzSingleLayerOperator>
          &incident_matrix,
      const std::complex<double> wavenumber, const int dofs_vector,
      const int dofs_scalar) {
    system_matrix->block(dofs_vector, 0, dofs_scalar, dofs_vector) =
        *P_times_mass_inverse *
        incident_matrix.get_incidence_matrix().transpose();
    return;
  }
  /**
   *  \brief Assemble the excitation block right of and below the M block.
   */
  static void assembleV(
      Eigen::MatrixXcd *system_matrix, Eigen::VectorXcd *excitation,
      const AnsatzSpace<HelmholtzSingleLayerOperator> &ansatz_space_scalar,
      const std::unique_ptr<Eigen::MatrixXcd> &P_times_mass_inverse,
      VoltageSource source, const int dofs_vector, const int dofs_scalar) {
    Eigen::MatrixXcd V = get_voltage_incidence_matrix(
        excitation, source, ansatz_space_scalar.get_number_of_dofs(),
        ansatz_space_scalar.get_polynomial_degree(),
        ansatz_space_scalar.get_refinement_level());

    system_matrix->block(dofs_vector, dofs_vector + dofs_scalar, dofs_scalar,
                         V.rows()) = *P_times_mass_inverse * V.transpose();
    system_matrix->block(dofs_vector + dofs_scalar, dofs_vector, V.rows(),
                         dofs_scalar) = V;
    return;
  }
  /**
   *  \brief Assemble the excitation block on the very right of and very bottom.
   *  Task of this block is to bundle all the degrees of freedom of one pole to
   *  the same potential.
   *
   *  This matrix has two columns each containing ones for the degrees of
   *  freedom corresponding either to the positive or negative pole of the
   *  voltage source. The first line mimes the difference of the potential.
   */
  static void assembleA(Eigen::MatrixXcd *system_matrix,
                        const VoltageSource &source, const int refinement_level,
                        const int dofs_vector, const int dofs_scalar) {
    const int dofs_excitation_positive =
        source.positive_.size() * (1 << refinement_level);
    const int dofs_excitation_negative =
        source.negative_.size() * (1 << refinement_level);
    const int dofs_excitation =
        dofs_excitation_positive + dofs_excitation_negative + 1;
    Eigen::MatrixXcd A = Eigen::MatrixXcd::Ones(dofs_excitation, 2);
    A.block(1, 1, dofs_excitation_negative, 1) =
        Eigen::VectorXcd::Zero(dofs_excitation_negative);
    A.block(1 + dofs_excitation_negative, 0, dofs_excitation_positive, 1) =
        Eigen::VectorXcd::Zero(dofs_excitation_positive);
    A(0, 0) = 1;
    A(0, 1) = -1;

    const int dofs_basis = dofs_vector + dofs_scalar;
    system_matrix->block(dofs_basis, dofs_basis + dofs_excitation,
                         dofs_excitation, 2) = A;
    system_matrix->block(dofs_basis + dofs_excitation, dofs_basis, 2,
                         dofs_excitation) = A.transpose();
    return;
  }
  static Eigen::MatrixXcd get_voltage_incidence_matrix(
      Eigen::VectorXcd *excitation, VoltageSource source, const int dofs_scalar,
      const int polynomial_degree, const int refinement_level) {
    const int elements_per_edge = (1 << refinement_level);
    const int num_elements_per_patch = elements_per_edge * elements_per_edge;
    const int dofs_excitation_V =
        (source.positive_.size() + source.negative_.size()) *
            elements_per_edge +
        1;
    const int number_of_basisfunctions = polynomial_degree + elements_per_edge;

    Eigen::MatrixXcd ret =
        Eigen::MatrixXcd::Zero(dofs_excitation_V, dofs_scalar);

    const int number_of_patches_negative = source.negative_.size();
    int row_count = 1;
    for (auto i = 0; i < number_of_patches_negative; ++i) {
      std::vector<int> edge_dofs = GlueRoutines::getEdgeDofIndices(
          source.negative_[i][1], number_of_basisfunctions,
          number_of_basisfunctions,
          source.negative_[i][0] * num_elements_per_patch);
      for (auto dof : edge_dofs) {
        ret(row_count, dof) = -1;
        excitation[0](dof) = -1;
        ++row_count;
      }
    }

    const int number_of_patches_positive = source.positive_.size();
    for (auto i = 0; i < number_of_patches_positive; ++i) {
      std::vector<int> edge_dofs = GlueRoutines::getEdgeDofIndices(
          source.positive_[i][1], number_of_basisfunctions,
          number_of_basisfunctions,
          source.positive_[i][0] * num_elements_per_patch);
      for (auto dof : edge_dofs) {
        ret(row_count, dof) = -1;
        excitation[0](dof) = 1;
        ++row_count;
      }
    }
    return ret;
  }
};
}  // namespace Bembel
#endif  // BEMBEL_SRC_AUGMENTEDEFIE_AUGMENTEDEFIEASSEMBLER_HPP_
