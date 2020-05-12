// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_ANSATZSPACE_PROJECTOR_H_
#define BEMBEL_ANSATZSPACE_PROJECTOR_H_

namespace Bembel {
namespace ProjectorRoutines {
template <typename Derived, unsigned int DF>
struct projector_matrixmaker_ {
  static Eigen::SparseMatrix<double> makeMatrix(
      const Bembel::SuperSpace<Derived>& super_space,
      const int knotrepetition) {
    assert(false && "This needs to be specialized");
    return Eigen::SparseMatrix<double>(1, 1);
  }
};
}  // namespace ProjectorRoutines

/**
 *  \ingroup AnsatzSpace
 *  \brief The projector solves interpolation problems to identify the the
 * correct linear combination to represent a B-Spline basis local to each patch
 * by means of the given superspace, and then assembles the required
 * transformation matrices.
 */
template <typename Derived>
class Projector {
 public:
  //////////////////////////////////////////////////////////////////////////////
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  Projector() {}
  Projector(const SuperSpace<Derived>& super_space, const int knot_repetition) {
    init_Projector(super_space, knot_repetition);
  }
  //////////////////////////////////////////////////////////////////////////////
  /// init
  //////////////////////////////////////////////////////////////////////////////
  void init_Projector(const SuperSpace<Derived>& super_space,
                      const int knot_repetition) {
    knot_repetition_ = knot_repetition;
    projector_ = makeProjectionMatrix(super_space, knot_repetition);
    dofs_before_projector_ = projector_.rows();
    dofs_after_projector_ = projector_.cols();
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  /// getter
  //////////////////////////////////////////////////////////////////////////////
  int get_knot_repetition() const { return knot_repetition_; }
  const Eigen::SparseMatrix<double>& get_projection_matrix() {
    return projector_;
  }
  int get_dofs_after_projector() const { return dofs_after_projector_; }
  //////////////////////////////////////////////////////////////////////////////
  /// private members
  //////////////////////////////////////////////////////////////////////////////
 private:
  Eigen::SparseMatrix<double> makeProjectionMatrix(
      const SuperSpace<Derived>& super_space, const int knotrepetition) {
    return ProjectorRoutines::projector_matrixmaker_<
        Derived,
        LinearOperatorTraits<Derived>::Form>::makeMatrix(super_space,
                                                         knotrepetition);
  }
  // we declare functionality which has not been implemented (yet)
  // to be private
  int knot_repetition_;
  int dofs_before_projector_;
  int dofs_after_projector_;
  Eigen::SparseMatrix<double> projector_;
  Projector(const Projector<Derived>& other);
  Projector(Projector<Derived>&& other);
  Projector& operator=(const Projector<Derived>& other);
  Projector& operator=(Projector<Derived>&& other);
  //////////////////////////////////////////////////////////////////////////////
};

namespace ProjectorRoutines {

// Local helper struct
struct _proj_info {
  std::vector<double> vals;
  std::vector<int> rows;
  std::vector<int> cols;
  int num_c_dof;
  int num_dc_dof;
};

template <typename Derived>
inline _proj_info makeLocalProjectionTriplets(
    const SuperSpace<Derived>& super_space, const int pp1x, const int pp1y,
    const int knotrepetition_in) {
  using namespace Spl;
  const int M = super_space.get_refinement_level();
  const int maximal_polynomial_degree = std::max(pp1x, pp1y);
  assert(maximal_polynomial_degree == super_space.get_polynomial_degree() + 1 &&
         "superSpace not suitable for desired ansatz space -> polynomial "
         "degree mismatch");
  const int minp = std::min(pp1x, pp1y);
  const int knotrepetition = (maximal_polynomial_degree <= knotrepetition_in)
                                 ? minp
                                 : knotrepetition_in;
  const int n = (1 << M);
  const int patch_number = super_space.get_number_of_patches();

  // At the moment, we allow a difference of one only, i.e., we start our de
  // Rham sequence with a space of the same degree in every TP-direction.
  assert(std::abs(maximal_polynomial_degree - minp) <= 1);

  // Now, we construct the continuous spaces which will be the preimage of the
  // projector. n-1 is passed, since n = number_elements but n-1 = number of
  // interior knots.
  std::vector<double> c_space_knot_x =
      MakeUniformKnotVector(pp1x, n - 1, knotrepetition);
  std::vector<double> c_space_knot_y =
      MakeUniformKnotVector(pp1y, n - 1, knotrepetition);
  const int c_space_dim_x = c_space_knot_x.size() - pp1x;
  const int c_space_dim_y = c_space_knot_y.size() - pp1y;
  const int c_space_dim = c_space_dim_x * c_space_dim_y;

  _proj_info out;

  out.cols.reserve(c_space_dim * maximal_polynomial_degree *
                   maximal_polynomial_degree * patch_number);
  out.vals.reserve(c_space_dim * maximal_polynomial_degree *
                   maximal_polynomial_degree * patch_number);
  out.rows.reserve(c_space_dim * maximal_polynomial_degree *
                   maximal_polynomial_degree * patch_number);

  std::vector<double> mask = MakeInterpolationMask(maximal_polynomial_degree);
  const int masksize = mask.size();

  // Now we assemble our system which we use for the interpolation
  assert(masksize == maximal_polynomial_degree &&
         "projector.cpp: System needs to be square");
  Eigen::Matrix<double, -1, -1> system(
      masksize * masksize,
      maximal_polynomial_degree * maximal_polynomial_degree);

  // Here, we suddenly use the degree for basisevaluation, i.e.,
  // maximal_polynomial_degree-1. This is confusing, but correct and tested.
  {
    double vals_y[Constants::MaxP + 1];
    double vals_x[Constants::MaxP + 1];
    for (int iy = 0; iy < maximal_polynomial_degree; ++iy) {
      Bembel::Basis::ShapeFunctionHandler::evalBasis(
          maximal_polynomial_degree - 1, vals_y, mask[iy]);
      for (int ix = 0; ix < maximal_polynomial_degree; ++ix) {
        Bembel::Basis::ShapeFunctionHandler::evalBasis(
            maximal_polynomial_degree - 1, vals_x, mask[ix]);
        for (int jy = 0; jy < maximal_polynomial_degree; ++jy) {
          for (int jx = 0; jx < maximal_polynomial_degree; ++jx) {
            system(iy * masksize + ix, jy * maximal_polynomial_degree + jx) =
                vals_x[jx] * vals_y[jy];
          }
        }
      }
    }
  }

  // Do LU-Decomposition of
  // systemsuper_space.get_mesh().get_element_tree().begin()
  const Eigen::PartialPivLU<Eigen::Matrix<double, -1, -1>> pplu(system);

  unsigned int element_number = 0;

  for (auto element = super_space.get_mesh().get_element_tree().cpbegin();
       element != super_space.get_mesh().get_element_tree().cpend();
       ++element) {
    // s = x
    const double pos_x = element->llc_(0);
    const double pos_y = element->llc_(1);

    const double h = element->get_h();
    const double mid_x = pos_x + (.5 * h);
    const double mid_y = pos_y + (.5 * h);

    // The following block of code computes the indeces of all continuous basis
    // functions which have a support on the element of index 'element'
    std::vector<double> nonzero_dofs_x;
    std::vector<double> nonzero_dofs_y;
    nonzero_dofs_x.reserve(pp1x);
    nonzero_dofs_y.reserve(pp1y);
    for (int dof_y = 0; dof_y < c_space_dim_y; ++dof_y) {
      if ((c_space_knot_y[dof_y] <= mid_y) &&
          (c_space_knot_y[dof_y + pp1y] >= mid_y)) {
        nonzero_dofs_y.push_back(dof_y);
      }
    }
    for (int dof_x = 0; dof_x < c_space_dim_x; ++dof_x) {
      if ((c_space_knot_x[dof_x] <= mid_x) &&
          (c_space_knot_x[dof_x + pp1x] >= mid_x)) {
        nonzero_dofs_x.push_back(dof_x);
      }
    }
    nonzero_dofs_x.shrink_to_fit();
    nonzero_dofs_y.shrink_to_fit();

    // We now loop over all basis functions wich are non-zero on the element
    for (auto dof_y : nonzero_dofs_y) {
      for (auto dof_x : nonzero_dofs_x) {
        // The next few lines are to check if the element is within the support
        // of the basis function given by dof_x and dof_y.
        std::vector<double> local_mask_x(masksize);
        std::vector<double> local_mask_y(masksize);
        for (int i = 0; i < masksize; ++i) {
          local_mask_x[i] = pos_x + (1.0 / n) * mask[i];
          local_mask_y[i] = pos_y + (1.0 / n) * mask[i];
          assert(local_mask_x[i] < 1 && local_mask_x[i] > 0);
          assert(local_mask_y[i] < 1 && local_mask_y[i] > 0);
        }

        // build unit vectors
        std::vector<double> c_coefs_x(c_space_dim_x, 0.0);
        std::vector<double> c_coefs_y(c_space_dim_y, 0.0);
        c_coefs_x[dof_x] = 1.;
        c_coefs_y[dof_y] = 1.;

        // evaluate rhs for interpolation problem
        std::vector<double> vals_x =
            DeBoor(c_coefs_x, c_space_knot_x, local_mask_x);
        std::vector<double> vals_y =
            DeBoor(c_coefs_y, c_space_knot_y, local_mask_y);

        // assemble and solve interpolation problem
        Eigen::VectorXd rhs(masksize * masksize);
        for (int iy = 0; iy < masksize; ++iy)
          for (int ix = 0; ix < masksize; ++ix)
            rhs(iy * masksize + ix) = vals_x[ix] * vals_y[iy];
        Eigen::VectorXd sol = pplu.solve(rhs);

        // insert entries into helperstruct
        for (int i = 0;
             i < maximal_polynomial_degree * maximal_polynomial_degree; ++i) {
          if (std::abs(sol(i)) > Constants::projector_tolerance) {
            out.cols.push_back(element->patch_ * (c_space_dim) +
                               (c_space_dim_x * dof_y) + dof_x);
            out.rows.push_back((element_number * maximal_polynomial_degree *
                                maximal_polynomial_degree) +
                               i);
            out.vals.push_back(sol(i));
          }
        }
      }
    }
    ++element_number;
  }

  out.num_c_dof = c_space_dim_x * c_space_dim_y * patch_number;
  out.num_dc_dof = maximal_polynomial_degree * maximal_polynomial_degree * n *
                   n * patch_number;

  out.cols.shrink_to_fit();
  out.rows.shrink_to_fit();
  out.vals.shrink_to_fit();

  assert((out.vals.size() == out.rows.size()) &&
         (out.rows.size() == out.cols.size()) &&
         "projector.cpp: you made a mistake, try again.");

  return out;
}

template <typename Derived>
struct projector_matrixmaker_<Derived, DifferentialForm::Continuous> {
  static Eigen::SparseMatrix<double> makeMatrix(
      const SuperSpace<Derived>& super_space, const int knotrepetition) {
    const int P = super_space.get_polynomial_degree();
    assert(P > 0 && "P must be 1 or larger for this type of discrete space");
    assert(knotrepetition <= P &&
           "Knot repetition must be smaller than P for this type of discrete "
           "space");
    // The matrices coincide with the 2 case, up to the limitations above
    Eigen::SparseMatrix<double> local_matrix = projector_matrixmaker_<
        Derived, DifferentialForm::Discontinuous>::makeMatrix(super_space,
                                                              knotrepetition);

    return local_matrix;
  }
};

template <typename Derived>
struct projector_matrixmaker_<Derived, DifferentialForm::DivConforming> {
  static Eigen::SparseMatrix<double> makeMatrix(
      const SuperSpace<Derived>& super_space, const int knotrepetition) {
    assert(super_space.get_polynomial_degree() > 0 &&
           "P must be 1 or larger for this type of discrete space");
    assert(knotrepetition <= super_space.get_polynomial_degree() &&
           "Knot repetition must be smaller than P for this type of discrete "
           "space");
    const auto info1 = ProjectorRoutines::makeLocalProjectionTriplets<Derived>(
        super_space, super_space.get_polynomial_degree() + 1,
        super_space.get_polynomial_degree(), knotrepetition);
    const auto info2 = ProjectorRoutines::makeLocalProjectionTriplets<Derived>(
        super_space, super_space.get_polynomial_degree(),
        super_space.get_polynomial_degree() + 1, knotrepetition);
    const int size1 = info1.vals.size();
    const int size2 = info2.vals.size();
    assert(size1 == size2);
    std::vector<Eigen::Triplet<double>> trips;
    for (int k = 0; k < size1; ++k) {
      trips.push_back(
          Eigen::Triplet<double>(info1.rows[k], info1.cols[k], info1.vals[k]));
    }
    for (int k = 0; k < size2; ++k) {
      trips.push_back(Eigen::Triplet<double>(info2.rows[k] + info1.num_dc_dof,
                                             info2.cols[k] + info1.num_c_dof,
                                             info2.vals[k]));
    }
    Eigen::SparseMatrix<double> local_matrix(
        info1.num_dc_dof + info2.num_dc_dof, info1.num_c_dof + info2.num_c_dof);
    local_matrix.setFromTriplets(trips.begin(), trips.end());

    return local_matrix;
  }
};

template <typename Derived>
struct projector_matrixmaker_<Derived, DifferentialForm::Discontinuous> {
  static Eigen::SparseMatrix<double> makeMatrix(
      const SuperSpace<Derived>& super_space, const int knotrepetition) {
    const auto info = ProjectorRoutines::makeLocalProjectionTriplets<Derived>(
        super_space, super_space.get_polynomial_degree() + 1,
        super_space.get_polynomial_degree() + 1, knotrepetition);
    const int size = info.vals.size();
    std::vector<Eigen::Triplet<double>> trips;
    for (int k = 0; k < size; ++k) {
      trips.push_back(
          Eigen::Triplet<double>(info.rows[k], info.cols[k], info.vals[k]));
    }
    Eigen::SparseMatrix<double> local_matrix(info.num_dc_dof, info.num_c_dof);
    local_matrix.setFromTriplets(trips.begin(), trips.end());

    return local_matrix;
  }
};
}  // namespace ProjectorRoutines

}  // namespace Bembel

#endif
