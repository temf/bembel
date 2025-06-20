// This file is part of Bembel, the higher order C++ boundary element library.
//
// Copyright (C) 2022 see <http://www.bembel.eu>
//
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_SRC_AUGMENTEDEFIE_AUGMENTEDEFIEEXCITATION_HPP_
#define BEMBEL_SRC_AUGMENTEDEFIE_AUGMENTEDEFIEEXCITATION_HPP_

namespace Bembel {
/**
 *  \ingroup AugmentedEFIE
 *  \brief This class, given a LinearForm with defined traits, takes care of the
 * assembly of the right hand side.
 */
template <typename Derived, typename LinOp>
class AugmentedEFIEExcitation {
 public:
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  AugmentedEFIEExcitation() {}
  explicit AugmentedEFIEExcitation(const AnsatzSpace<LinOp> &ansatz_space,
                                   const int dofs_scalar) {
    init_AugmentedEFIEExcitation(ansatz_space, dofs_scalar);
  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_AugmentedEFIEExcitation
  //////////////////////////////////////////////////////////////////////////////
  void init_AugmentedEFIEExcitation(const AnsatzSpace<LinOp> &ansatz_space,
                                    const int dofs_scalar) {
    ansatz_space_ = ansatz_space;
    disc_lf_ = DiscreteLinearForm<Derived, LinOp>(ansatz_space_);
    dofs_scalar_ = dofs_scalar;
    dofs_vector_ = ansatz_space.get_number_of_dofs();
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    compute
  //////////////////////////////////////////////////////////////////////////////
  /**
   *    \todo Add inline commentary
   **/
  void compute() {
    excitation_ = Eigen::VectorXcd::Zero(dofs_scalar_ + dofs_vector_);
    disc_lf_.compute();

    excitation_.head(dofs_vector_) = -disc_lf_.get_discrete_linear_form();
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    getter
  //////////////////////////////////////////////////////////////////////////////
  Derived &get_linear_form() { return disc_lf_.get_linear_form(); }
  const Eigen::VectorXcd &get_excitation() const { return excitation_; }
  //////////////////////////////////////////////////////////////////////////////
  //    private member variables
  //////////////////////////////////////////////////////////////////////////////
 private:
  Eigen::VectorXcd excitation_;

  DiscreteLinearForm<Derived, LinOp> disc_lf_;
  AnsatzSpace<LinOp> ansatz_space_;

  int dofs_scalar_;
  int dofs_vector_;
};

}  // namespace Bembel
#endif  // BEMBEL_SRC_AUGMENTEDEFIE_AUGMENTEDEFIEEXCITATION_HPP_
