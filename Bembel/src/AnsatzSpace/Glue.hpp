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
#ifndef BEMBEL_SRC_ANSATZSPACE_GLUE_HPP_
#define BEMBEL_SRC_ANSATZSPACE_GLUE_HPP_

namespace Bembel {
/*

Every edge has an index between 0 and 3, where
0 = (0,0)->(1,0)
1 = (1,0)->(1,1)
2 = (1,1)->(0,1)
3 = (0,1)->(0,0)
*/

/*
  Example of Dof enumeration on one patch:
  dimXdir = 4;
  dimYdir = 3;
                  Edge 2      (1,1)

    y         |8  9  10 11|
    A  Edge 3 |4  5  6  7 |  Edge 1
    |         |0  1  2  3 |

                  Edge 0
   (0,0)  -> x
*/

namespace GlueRoutines {
/**
 * \brief Helper struct.
 */
struct dofIdentification {
  std::vector<int> dofs;
  int coef;
};
/**
 * \brief Helper struct for assembling the Glue for the different cases.
 */
template <typename Derived, unsigned int DF>
struct glue_identificationmaker_ {
  static std::vector<dofIdentification> makeIdentification(
      const std::vector<std::array<int, 4>> &edges_,
      const SuperSpace<Derived> &superspace, const Projector<Derived> &proj) {
    assert(false && "Needs to be specialized");
    return {};
  }
};

}  // namespace GlueRoutines

/**
 * \ingroup AnsatzSpace
 * \brief This class takes care of identifying DOFs on different edges, which
 *must be identified with one another.
 *
 * The Glue class has routines to assemble continuous B-splines across patch
 *boundaries depending on the template parameter
 *LinearOperatorTraits<Derived>::Form. If this template parameter is equal to
 *DifferentialForm::Continuous, the B-splines are globally glued continuously.
 *Otherwise, it can also be DifferentialForm::DivConforming, in which case the
 *normal component is glued continuously. Finally, there is the option of not
 *gluing the B-splines at all if the parameter is
 *DifferentialForm::Discontinuous.
 **/
template <typename Derived>
class Glue {
  std::vector<std::array<int, 4>> edges_;
  int dofs_after_glue;
  Eigen::SparseMatrix<double> glue_mat_;

 public:
  /**
   * \brief Constructor for the Glue class.
   *
   * This constructor initializes a Glue object with the provided SuperSpace and
   * Projector objects.
   *
   * \param superspace The SuperSpace to handle the basis functions.
   * \param proj The Projector which provides information on the dofs before
   *        gluing.
   */
  Glue(const SuperSpace<Derived> &superspace, const Projector<Derived> &proj) {
    init_Glue(superspace, proj);
  }
  /**
   * \brief Initializes the AnsatzSpace object.
   *
   * This function initializes the member variables of the Glue object.
   *
   * \param superspace The SuperSpace to handle the basis functions.
   * \param proj The Projector which provides information on the dofs before
   *        gluing.
   */
  void init_Glue(const SuperSpace<Derived> &superspace,
                 const Projector<Derived> &proj) {
    edges_ = superspace.get_mesh().get_element_tree().patchTopologyInfo();
    glue_mat_ = assembleGlueMatrix(superspace, proj);
    return;
  }
  /**
   * \brief Returns Glue matrix to assemble continuous B-splines.
   *
   * This function returns the Glue matrix which transforms smooth B-splines to
   * continuous shape functions across patch boundaries. If continuity is not
   * needed this matrix is the identity. Other wise this matrix is a tall thin
   * matrix combining degrees of freedom.
   *
   * \return Eigen::SparseMatrix<double> Glue matrix
   */
  Eigen::SparseMatrix<double> get_glue_matrix() const { return glue_mat_; }
  /**
   * \brief Generates a list of degrees of freedom (DOFs) identifications.
   *
   * This function creates a list of degrees of freedom (DOFs) identifications
   * which need to be glued together. The return type is a std::vector of
   * GlueRoutines::dofIdentification which collects the DOFs to be glued.
   * Furthermore, the coef in this struct denotes the orientation.
   *
   * \param superspace The SuperSpace reference to handle basis functions.
   * \param proj The Projector contains the number of DOFs.
   * \return A vector of GlueRoutines::dofIdentification containing the DOF
   * identifications.
   */
  inline std::vector<GlueRoutines::dofIdentification> makeDofIdentificationList(
      const SuperSpace<Derived> &superspace, const Projector<Derived> &proj) {
    return GlueRoutines::glue_identificationmaker_<
        Derived,
        LinearOperatorTraits<Derived>::Form>::makeIdentification(edges_,
                                                                 superspace,
                                                                 proj);
  }
  /**
   * \brief Assembles the Glue matrix according to the DOF identification list.
   *
   * This function first sorts the DOFs in each identification list. Than the
   * identification lists get sorted according to the first entries. After that
   * the DOFs get stored in the Glue matrix.
   *
   * \param superspace The SuperSpace reference to handle basis functions.
   * \param proj The Projector contains the number of DOFs.
   * \return The Glue matrix.
   */
  Eigen::SparseMatrix<double> assembleGlueMatrix(
      const SuperSpace<Derived> &superspace, const Projector<Derived> &proj) {
    const int pre_dofs = proj.get_dofs_after_projector();
    // The dofs that need to be identified with each other are divided into
    // master and slaves, where the master is the firs w.r.t. the tp-ordering
    std::vector<Eigen::Triplet<double>> trips;
    std::vector<GlueRoutines::dofIdentification> dof_id =
        makeDofIdentificationList(superspace, proj);

    // We keep track of certain information
    std::vector<bool> dof_is_slave(pre_dofs, false);
    std::vector<bool> dof_is_master(pre_dofs, false);

    // Here, we fill the two vectors above
    int number_of_slaves = 0;
    for (auto dofset : dof_id) {
      // This sorting is required, since by construction only dofs[0]<dofs[1] is
      // given, but in the case of corners and H1 discretizations, it might
      // happen that dofs[0]>dofs[2]. Moreover, we count how many dofs we "glue
      // away".
      std::sort(dofset.dofs.begin(), dofset.dofs.end());
      dof_is_master[dofset.dofs[0]] = true;
      for (int i = 1; i < dofset.dofs.size(); ++i) {
        dof_is_slave[dofset.dofs[i]] = true;
        ++number_of_slaves;
      }
    }

    // Now, we sort w.r.t. the masters.
    std::sort(dof_id.begin(), dof_id.end(),
              [](GlueRoutines::dofIdentification a,
                 GlueRoutines::dofIdentification b) {
                return a.dofs[0] < b.dofs[0];
              });

    const int post_dofs = pre_dofs - number_of_slaves;

    // This block is the heart of the algorithms. Skip keeps track of how many
    // slaves have been skipped already, and master_index keeps track on which
    // master is about to be glued next.
    int skip = 0;
    int master_index = 0;
    for (int i = 0; i < post_dofs; ++i) {
      const int post_index = i;
      while (dof_is_slave[i + skip] && ((i + skip) < pre_dofs)) ++skip;
      const int pre_index = i + skip;
      trips.push_back(Eigen::Triplet<double>(pre_index, post_index, 1));
      // The dof cannot be declared slave and master at the same time
      assert(!(dof_is_master[pre_index] && dof_is_slave[pre_index]));
      if (dof_is_master[pre_index]) {
        // The dofs in dof_id don't know they might be moved forward. So the
        // smallest dof of those to be identified should coincide with the
        // pre_index.
        assert(pre_index == dof_id[master_index].dofs[0]);
        const int number_of_partners = dof_id[master_index].dofs.size();
        for (int j = 1; j < number_of_partners; ++j) {
          trips.push_back(Eigen::Triplet<double>(dof_id[master_index].dofs[j],
                                                 post_index,
                                                 dof_id[master_index].coef));
        }
        ++master_index;
      }
    }

    Eigen::SparseMatrix<double> glue_matrix(pre_dofs, post_dofs);
    glue_matrix.setFromTriplets(trips.begin(), trips.end());

    assert(trips.size() == pre_dofs);
    return glue_matrix;
  }
};

namespace GlueRoutines {

/**
 * \brief This routine collects the indices of DOFs with support on an edge.
 *
 * The indices are computed depending on the edge case and number of DOFs in the
 * x and y direction of the tensor product. The shift parameter can be used to
 * shift all dofs according to the patch index or vector component.
 *
 * Edge Case
 *
 * 0 = (0,0)->(1,0)
 *
 * 1 = (1,0)->(1,1)
 *
 * 2 = (1,1)->(0,1)
 *
 * 3 = (0,1)->(0,0)
 *
 * \param edgeCase Edge case of the patch to collect the DOFs.
 * \param dimXdir Number of DOFs in the x direction of the tensor product.
 * \param dimYdir Number of DOFs in the y direction of the tensor product.
 * \param shift Offset of the DOFs.
 */
std::vector<int> getEdgeDofIndices(int edgeCase, int dimXdir, int dimYdir,
                                   int shift) {
  std::vector<int> out;
  switch (edgeCase) {
    case (0): {
      out.reserve(dimXdir);
      assert(
          "This should be a given. If not something went wrong in Glue.hpp" &&
          dimXdir <= dimYdir);
      for (int i = 0; i < dimXdir; ++i) {
        out.push_back(i + shift);
      }
      return out;
    };
    case (1): {
      out.reserve(dimYdir);
      assert(
          "This should be a given. If not something went wrong in Glue.hpp" &&
          dimYdir <= dimXdir);
      for (int i = 0; i < dimYdir; ++i) {
        out.push_back(dimXdir * (i + 1) - 1 + shift);
      }
      return out;
    };
    case (2): {
      out.reserve(dimXdir);
      assert(
          "This should be a given. If not something went wrong in Glue.hpp" &&
          dimXdir <= dimYdir);
      for (int i = 0; i < dimXdir; ++i) {
        out.push_back(dimXdir * (dimYdir - 1) + i + shift);
      }
      return out;
    };
    case (3): {
      out.reserve(dimYdir);
      assert(
          "This should be a given. If not something went wrong in Glue.hpp" &&
          dimYdir <= dimXdir);
      for (int i = 0; i < dimYdir; ++i) {
        out.push_back(dimXdir * i + shift);
      }
      return out;
    };
    default: {
    };
      // An edge might have a -1 index. This occurs only when no partner could
      // be found, and the -1 is the placeholder of the missing partner.
  }
  return out;
}

// The following 7 functions figure out how edges and the vector components
// are oriented w.r.t. each other.
inline bool edgeIsForwardParametrized(int edgeCase) {
  return edgeCase == 0 || edgeCase == 1;
}
inline bool edgeIsBackwardsParametrized(int edgeCase) {
  return !(edgeIsForwardParametrized(edgeCase));
}
inline bool normalComponentIsInwardDirected(int edgeCase) {
  return edgeCase == 0 || edgeCase == 3;
}
inline bool normalComponentIsOutwardDirected(int edgeCase) {
  return !(normalComponentIsInwardDirected(edgeCase));
}

inline bool reverseParametrized(const std::array<int, 4> &edge) {
  return ((edgeIsForwardParametrized(edge[2]) &&
           edgeIsForwardParametrized(edge[3])) ||
          (edgeIsBackwardsParametrized(edge[2]) &&
           edgeIsBackwardsParametrized(edge[3])))
             ? true
             : false;
}

inline int glueCoefficientDivergenceConforming(const std::array<int, 4> &edge) {
  return ((normalComponentIsInwardDirected(edge[2]) &&
           normalComponentIsInwardDirected(edge[3])) ||
          (normalComponentIsOutwardDirected(edge[2]) &&
           normalComponentIsOutwardDirected(edge[3])))
             ? -1
             : 1;
}
inline bool edgeToBeGluedInFirstComp(int edgeCase) {
  return (edgeCase == 0 || edgeCase == 2) ? false : true;
}

/**
 * \brief Specification of the struct for the discontinuous case.
 *
 * Nothing needs to be glued.
 */
template <typename Derived>
struct glue_identificationmaker_<Derived, DifferentialForm::Discontinuous> {
  static std::vector<dofIdentification> makeIdentification(
      const std::vector<std::array<int, 4>> &edges_,
      const SuperSpace<Derived> &superspace, const Projector<Derived> &proj) {
    return {};
  }
};

/**
 * \brief Specification of the struct for the scalar continuous case.
 *
 * The scalar case: Here, we have only one vector component to worry about.
 * However, we need to take care of edges, since here multiple (more than 2)
 * dofs will be reduced to one. This is also the reason we we work with the
 * dofIdentification struct and not already with Eigen::Triplet. The function
 * builds upon the fact, that px == py == p.
 */
template <typename Derived>
struct glue_identificationmaker_<Derived, DifferentialForm::Continuous> {
  static std::vector<dofIdentification> makeIdentification(
      const std::vector<std::array<int, 4>> &edges_,
      const SuperSpace<Derived> &superspace, const Projector<Derived> &proj) {
    // We check if the space can even be continuous globally.
    assert(superspace.get_polynomial_degree() >= 1);
    const int one_d_dim = superspace.get_polynomial_degree() + 1 +
                          proj.get_knot_repetition() *
                              ((1 << superspace.get_refinement_level()) - 1);
    const int dofs_per_patch = one_d_dim * one_d_dim;

    std::vector<GlueRoutines::dofIdentification> out;
    out.reserve(edges_.size() * one_d_dim);

    std::vector<int> already_stored_in(proj.get_dofs_after_projector(), -1);
    int d_count = 0;
    for (auto edge : edges_) {
      // The shift is to account for different patches, since the
      // getEdgeDof-routines only can enumerate w.r.t. patch 0.
      const int shift_e1 = edge[0] * dofs_per_patch;
      const int shift_e2 = edge[1] * dofs_per_patch;
      const bool needReversion = reverseParametrized(edge);

      std::vector<int> dofs_e1 =
          getEdgeDofIndices(edge[2], one_d_dim, one_d_dim, shift_e1);
      std::vector<int> dofs_e2 =
          getEdgeDofIndices(edge[3], one_d_dim, one_d_dim, shift_e2);

      // Check if the edge is hanging. For screens one would need an "else".
      if (edge[1] > -1 && edge[0] > -1) {
        for (int i = 0; i < one_d_dim; ++i) {
          GlueRoutines::dofIdentification d;
          const int j = needReversion ? one_d_dim - 1 - i : i;
          // const int j = i;
          assert(dofs_e1[i] != dofs_e2[j] &&
                 "If this happens something went horribly wrong.");

          const int small_dof = std::min(dofs_e1[i], dofs_e2[j]);
          const int large_dof = std::max(dofs_e1[i], dofs_e2[j]);

          // In the continuous case, more than two dofs might be glued together.
          // Therefore, we must check if we know one of the dofs already. If
          // yes, we just add the new one to the dof_id, if not, we make a new
          // identification group.
          if (already_stored_in[small_dof] > -1) {
            // It could be that the large_dof is already matched with another
            // dof, i.e., in the case of a 'circle'. Then there is nothing to
            // do. If not, we match it with the master of the partner.
            if (already_stored_in[large_dof] == -1) {
              out[already_stored_in[small_dof]].dofs.push_back(large_dof);
              already_stored_in[large_dof] = already_stored_in[small_dof];
            } else {
              if (!(already_stored_in[small_dof] ==
                    already_stored_in[large_dof])) {
                // This case is tricky. Assume that we have to identify four
                // DOFs with each other, but they have already been assigned in
                // pairs of two. Then we need to reverse this process. First, we
                // grab the two storage locations.
                const int small_store = std::min(already_stored_in[large_dof],
                                                 already_stored_in[small_dof]);
                const int large_store = std::max(already_stored_in[large_dof],
                                                 already_stored_in[small_dof]);
                for (auto dfs : out[large_store].dofs) {
                  // now we put all of those in the larger location into the
                  // smaller one.
                  already_stored_in[dfs] = small_store;
                  for (auto other_dfs : out[small_store].dofs) {
                    assert(dfs != other_dfs);
                  }
                  out[small_store].dofs.push_back(dfs);
                }
                // Now we set the larger storage location to empty.
                out[large_store].dofs = {};
                out[large_store].coef = 0;
              }
            }
          } else if (already_stored_in[large_dof] > -1) {
            // It could be that the small_dof is already matched with another
            // dof, i.e., in the case of a 'circle'. Then there is nothing to
            // do. If not, we match it with the master of the partner.
            if (already_stored_in[small_dof] == -1) {
              out[already_stored_in[large_dof]].dofs.push_back(small_dof);
              already_stored_in[small_dof] = already_stored_in[large_dof];
            } else {
              if (!(already_stored_in[small_dof] ==
                    already_stored_in[large_dof])) {
                // This case is tricky. Assume that we have to identify four
                // DOFs with each other, but they have already been assigned in
                // pairs of two. Then we need to reverse this process. First, we
                // grab the two storage locations.
                const int small_store = std::min(already_stored_in[large_dof],
                                                 already_stored_in[small_dof]);
                const int large_store = std::max(already_stored_in[large_dof],
                                                 already_stored_in[small_dof]);
                for (auto dfs : out[large_store].dofs) {
                  // now we put all of those in the larger location into the
                  // smaller one.
                  already_stored_in[dfs] = small_store;
                  for (auto other_dfs : out[small_store].dofs) {
                    assert(dfs != other_dfs);
                  }
                  out[small_store].dofs.push_back(dfs);
                }
                // Now we set the larger storage location to empty.
                out[large_store].dofs = {};
                out[large_store].coef = 0;
              }
            }
          } else {
            // With the exception of corners, this will be the default case.
            // We just add a pair of dofs to the bookkeeping.
            d.dofs.push_back(small_dof);
            d.dofs.push_back(large_dof);
            already_stored_in[small_dof] = d_count;
            already_stored_in[large_dof] = d_count;
            d_count++;
            d.coef = 1;
            out.push_back(d);
          }
        }
      }
    }
    // Now we need to clean up the empty dofsets, since subsequent routines
    // assume there to be at least one element.
    for (auto x = out.begin(); x != out.end(); ++x) {
      if ((*x).dofs.size() == 0) {
        out.erase(x);
      }
    }
    return out;
  }
};

/**
 * \brief Specification of the struct for the vector Div conforming case.
 *
 * The Maxwell case: Here, the identification is 1-to-1, but we need to take
 * care to glue the right vector component. The function builds upon the fact,
 * that px == py == p.
 */
template <typename Derived>
struct glue_identificationmaker_<Derived, DifferentialForm::DivConforming> {
  static std::vector<dofIdentification> makeIdentification(
      const std::vector<std::array<int, 4>> &edges_,
      const SuperSpace<Derived> &superspace, const Projector<Derived> &proj) {
    // since we assume px == py and uniform refinement in both directions, there
    // will be a small_dim and a large_dim in every vector component.
    const int small_dim = superspace.get_polynomial_degree() +
                          proj.get_knot_repetition() *
                              ((1 << superspace.get_refinement_level()) - 1);
    const int large_dim = superspace.get_polynomial_degree() + 1 +
                          proj.get_knot_repetition() *
                              ((1 << superspace.get_refinement_level()) - 1);
    const int dofs_per_patch_per_component = small_dim * large_dim;
    // sanity check
    assert(dofs_per_patch_per_component ==
               (proj.get_dofs_after_projector() /
                (superspace.get_number_of_patches() * 2)) &&
           "The assembly of the glue matrix is highly specific to the space. "
           "Something went wrong; is the discrete space correct?");

    std::vector<GlueRoutines::dofIdentification> out;
    out.reserve(edges_.size() * small_dim);

    for (auto edge : edges_) {
      assert(edge[0] <= edge[1] || edge[1] == -1);
      // The DOFs that need to be glued on edges 0 and 2 are in the second
      // vector component. The remainder of the dof-index-shift is to account
      // for dofs on other patches.
      const int shift_e1 = edge[0] * dofs_per_patch_per_component +
                           (edgeToBeGluedInFirstComp(edge[2])
                                ? 0
                                : (dofs_per_patch_per_component *
                                   superspace.get_number_of_patches()));
      const int shift_e2 = edge[1] * dofs_per_patch_per_component +
                           (edgeToBeGluedInFirstComp(edge[3])
                                ? 0
                                : (dofs_per_patch_per_component *
                                   superspace.get_number_of_patches()));

      const int xdir_dim_1 =
          edgeToBeGluedInFirstComp(edge[2]) ? large_dim : small_dim;
      const int ydir_dim_1 =
          edgeToBeGluedInFirstComp(edge[2]) ? small_dim : large_dim;
      const int xdir_dim_2 =
          edgeToBeGluedInFirstComp(edge[3]) ? large_dim : small_dim;
      const int ydir_dim_2 =
          edgeToBeGluedInFirstComp(edge[3]) ? small_dim : large_dim;
      std::vector<int> dofs_e1 =
          getEdgeDofIndices(edge[2], xdir_dim_1, ydir_dim_1, shift_e1);
      std::vector<int> dofs_e2 =
          getEdgeDofIndices(edge[3], xdir_dim_2, ydir_dim_2, shift_e2);

      const int size_of_edge_dofs = dofs_e1.size();
      assert(size_of_edge_dofs == dofs_e2.size() || edge[1] == -1);

      const bool needReversion = reverseParametrized(edge);
      const int coef = glueCoefficientDivergenceConforming(edge);

      // Check if the edge is hanging. For screens one would need an "else".
      if (edge[1] > -1 && edge[0] > -1) {
        for (int i = 0; i < size_of_edge_dofs; ++i) {
          GlueRoutines::dofIdentification d;
          const int j = needReversion ? size_of_edge_dofs - 1 - i : i;
          assert(dofs_e1[i] != dofs_e2[j] &&
                 "If this happens something went horribly wrong.");
          d.dofs.push_back(std::min(dofs_e1[i], dofs_e2[j]));
          d.dofs.push_back(std::max(dofs_e1[i], dofs_e2[j]));
          d.coef = coef;
          out.push_back(d);
        }
      }
    }
    out.shrink_to_fit();
    return out;
  }
};
};  // namespace GlueRoutines

}  // namespace Bembel
#endif  // BEMBEL_SRC_ANSATZSPACE_GLUE_HPP_
