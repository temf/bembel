// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "Data.hpp"
#include "bemtest.h"

using namespace Bembel;
/**
 *  @brief         This is a test to see if the HMatrix Wrappers and converters
 * are doing the same as the legacy stuff
 *
 */
int Test::test_HMatrixWrappert() {
  Geometry myGeom(Bembel::Test::mkScreen());
  MaxwellSingle myMax(std::complex<double>(1, 0));

  Discretization<MaxwellSingle> myDisc;

  const std::function<Eigen::Vector3cd(Eigen::Vector3d, std::complex<double>)>
      fun = [](Eigen::Vector3d pt, std::complex<double> kappa) {
        const Eigen::Vector3d position(0.2, 0.2, 0.2);
        const Eigen::Vector3d length(0, 0.1, 0.1);
        return Data::Dipole(pt, kappa, position, length);
      };

  const int np_max = 30;
  const double eta = 0.0016;

  for (auto M : {0})
    for (auto P : {1, 2, 3, 4}) {
      myDisc.init_Discretization(myGeom, myMax, P + 1, 1, M);

      // Eigen::VectorXcd rhs = get_rhs(myDisc, fun);
      // double * rhsptr = eigen2cmplxptr(rhs);

      auto hmatset = get_hmatrixsettings(np_max, &(myDisc.get_disc()));
      hmatrixfactory hmatfac;
      hmatset.eta = eta;
      hmatfac.disc = &(myDisc.get_disc());
      hmatfac.hmatset = &hmatset;

      ct_root *legacy_H;
      legacy_H = (ct_root *)calloc(16, sizeof(ct_root));

      init_cluster_tree(&hmatfac, &legacy_H);

      Eigen::HierarchicalMatrix<MaxwellSingle> myH(myDisc, np_max, eta);

      Eigen::VectorXcd rand = Eigen::VectorXcd::Random(myH.rows());
      double *randptr = eigen2cmplxptr(rand);
      double *outptr = (double *)calloc(sizeof(double), rand.rows() * 2);
      Eigen::VectorXcd out = myH * rand;
      H2l2_HtimesVsmallMaxwell(legacy_H, randptr, outptr);
      Eigen::VectorXcd legacy_out = cmplxptr2eigen(outptr, rand.rows() * 2);
      free(outptr);
      free(randptr);
      free(legacy_H);

      if ((legacy_out - out).norm() > 1e-14) {
        // std::cout << "Norm is : " << (legacy_out-out).norm() << "\n";
        return 1;
      }
    }

  return (0);
}
