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

#include <Bembel/AnsatzSpace>
#include <Bembel/Geometry>
#include <Bembel/IO>
#include <unsupported/Bembel/AugmentedEFIE>

#include "examples/Data.hpp"
#include "examples/Error.hpp"
#include "examples/Grids.hpp"

template <typename t>
class logspace {
 private:
  t curvalue, base;

 public:
  logspace(t first, t base) : curvalue(first), base(base) {}

  t operator()() {
    t retval = curvalue;
    curvalue *= base;
    return retval;
  }
};

int main() {
  using namespace Bembel;
  using namespace Eigen;

  const int refinement_level = 1;
  const int polynomial_degree = 1;

  const double f_min = 250e6;
  const double base = 1.0064;
  const int N = 1000;
  std::vector<double> frequencies;
  std::generate_n(std::back_inserter(frequencies), N,
                  logspace<double>(f_min, base));

  std::cout << "Number of Frequency Samples: " << frequencies.size()
            << std::endl;

  IO::Logger<10> logger("Balun_M" + std::to_string(refinement_level) + ".txt",
                        ",");
  logger.file("frequency", "Zmag", "Zphase");

  Geometry geometry("balun_spiral.dat");

  VoltageSource source;

  source.positive_ = {{4, 1}, {5, 1}, {6, 1}, {7, 1}, {8, 1}, {9, 1}};
  source.negative_ = {{0, 1}, {1, 1}, {2, 2}, {3, 2}};

  AnsatzSpace<InductanceMatrix> ansatz_space_vector(geometry, refinement_level,
                                                    polynomial_degree);

  AugmentedEFIE<Eigen::MatrixXcd, InductanceMatrix> augmented_efie(
      ansatz_space_vector, geometry);

  std::cout << "\n" << std::string(60, '=') << "\n";
  for (auto f = frequencies.begin(); f != frequencies.end(); ++f) {
    std::cout << std::string(60, '=') + "\nStarting Sample: "
              << std::to_string(f - frequencies.begin()) << " @ " << *f << "Hz"
              << std::endl;
    double omega = 2 * BEMBEL_PI * *f;
    std::complex<double> wavenumber = std::complex<double>(
        omega * std::sqrt(Constants::mu0 * Constants::eps0), 0.0);

    augmented_efie.set_wavenumber(wavenumber);
    augmented_efie.compute(source);

    const int dofs_vector = augmented_efie.get_dofs_vector();
    const int dofs_scalar = augmented_efie.get_dofs_scalar();
    const int dofs_total = augmented_efie.get_system_matrix().cols();

    // Define Excitation by Voltage Source
    Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(dofs_total);
    rhs(dofs_vector + dofs_scalar) = 1;

    // Solve System
    PartialPivLU<MatrixXcd> lu;
    lu.compute(augmented_efie.get_system_matrix());
    VectorXcd u = lu.solve(rhs);
    VectorXcd j = u.head(dofs_vector);

    int correction = (1 << refinement_level) * (1 << refinement_level);
    // Compute Impedance
    std::complex<double> Z = (u(dofs_total - 2) - u(dofs_total - 1)) /
                             u(dofs_vector + dofs_scalar) *
                             (double)(correction);

    std::cout << "Impedance: " << std::abs(Z) << std::endl;

    logger.file(*f, std::abs(Z), std::arg(Z) * 180 / BEMBEL_PI);
  }
  std::cout << "\n" << std::string(60, '=') << "\n";
  return 0;
}
