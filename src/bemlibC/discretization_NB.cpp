// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "discretization.h"

namespace Bembel {
discretization get_discretization_NB(const int deg, const int kntrep,
                                     pdeproblem *pde, meshdata *mesh) {
  discretization disc;

  /*
   * object stuff
   */
  disc.pde = pde;
  disc.mesh = mesh;

  // init_projectorr_0P(&disc);
  /*
   * fix constants depending on order of basis functions
   */
  disc.a_o = deg + 1;
  disc.a_bs = disc.a_o * disc.a_o;
  disc.a_bs2 = disc.a_bs * disc.a_bs;

  disc.na = init_projector(&disc, mesh->E.patch[0], mesh->M, disc.a_o, kntrep);
  // disc.na = mesh->nf;

  disc.real_na = pde->pdetype == Laplace ? disc.na : disc.na / 2;

  /*
   * quadrature related stuff
   */

  disc.quadrature_accuracy = pde->quadrature_accuracy(disc.a_o);

  disc.g_far = pde->g_far(disc.a_o);
  disc.g_pot = pde->g_pot(disc.a_o);

  /*
   * assign polynomials
   */

  disc.phi = select_phi(disc.a_o);                      //&phi_2B;
  disc.phi_dx = select_phi_dx(disc.a_o);                //&phi_2B;
  disc.phiphi = select_phiphi(disc.a_o);                //&phiphi_2B;
  disc.phiphi_dx = select_phiphi_dx(disc.a_o);          //&phiphi_dx_2B;
  disc.phiphi_dy = select_phiphi_dy(disc.a_o);          //&phiphi_dy_2B;
  disc.Phi_times_Phi = select_Phi_times_Phi(disc.a_o);  //&Phi_times_Phi_2B;
  if (pde->pdetype == Maxwell) {
    assert(disc.a_o > 1 && "Needs to be conforming");
  }
  disc.VPhi_scal_VPhi = select_VPhi_scal_bla(disc.a_o);
  disc.Div_Phi_times_Div_Phi = select_Div_Phi_times_Div_Phi(disc.a_o);
  disc.Curl_Phi_times_Curl_Phi = NULL;

#ifdef _BEMBEL_PRINT_INFO_
  printf("                           Ansatz functions of order %d\n", disc.a_o);
  printf("                           %d in total\n", disc.na);
#endif

  return disc;
}

discretization get_discretization_NB_Laplace_cont(const int deg,
                                                  const int kntrep,
                                                  pdeproblem *pde,
                                                  meshdata *mesh) {
  discretization disc;

  assert(deg > 1 && "Needs to be conforming");

  /*
   * object stuff
   */
  disc.pde = pde;
  disc.mesh = mesh;

  // init_projectorr_0P(&disc);
  /*
   * fix constants depending on order of basis functions
   */
  disc.a_o = deg;
  disc.a_bs = disc.a_o * disc.a_o;
  disc.a_bs2 = disc.a_bs * disc.a_bs;

  disc.na = init_projector_Laplace_cont(&disc, mesh->E.patch[0], mesh->M,
                                        disc.a_o, kntrep);
  // disc.na = mesh->nf;

  /*
   * quadrature related stuff
   */
  disc.quadrature_accuracy = pde->quadrature_accuracy(disc.a_o);
  disc.g_far = pde->g_far(disc.a_o);
  disc.g_pot = pde->g_pot(disc.a_o);

  /*
   * assign polynomials
   */

  disc.phi = select_phi(disc.a_o);                      //&phi_2B;
  disc.phi_dx = select_phi_dx(disc.a_o);                //&phi_2B;
  disc.phiphi = select_phiphi(disc.a_o);                //&phiphi_2B;
  disc.phiphi_dx = select_phiphi_dx(disc.a_o);          //&phiphi_dx_2B;
  disc.phiphi_dy = select_phiphi_dy(disc.a_o);          //&phiphi_dy_2B;
  disc.Phi_times_Phi = select_Phi_times_Phi(disc.a_o);  //&Phi_times_Phi_2B;
  disc.VPhi_scal_VPhi = select_VPhi_scal_bla(disc.a_o);
  disc.Div_Phi_times_Div_Phi = select_Div_Phi_times_Div_Phi(disc.a_o);
  disc.Curl_Phi_times_Curl_Phi = NULL;

#ifdef _BEMBEL_PRINT_INFO_
  printf("                           Ansatz functions of order %d\n", disc.a_o);
  printf("                           %d in total\n", disc.na);
#endif

  return disc;
}
}  // namespace Bembel