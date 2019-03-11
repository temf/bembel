// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "spline/patch.h"

namespace Spl {

void Patch::initHom(const std::vector<Eigen::Matrix<double, -1, -1>> &mats,
                    const std::vector<double> &knotx,
                    const std::vector<double> &knoty) {
  assert(mats.size() == 4);
  const int sxmat = mats[0].cols();
  const int symat = mats[0].rows();
  xuniqueknt = extract_unique_knots(knotx);
  yuniqueknt = extract_unique_knots(knoty);
  const int xnumpatch = xuniqueknt.size() - 1;
  const int ynumpatch = yuniqueknt.size() - 1;
  xp = knotx.size() - sxmat;
  yp = knoty.size() - symat;
  data.resize(4 * (xp * xnumpatch * yp * ynumpatch));
  {
    // Since its only for initialization, I do not care about speed.
    // Here I look weather the data given is already in bezier form.
    if (xuniqueknt.size() == 2 && yuniqueknt.size() == 2) {
      for (int i = 0; i < 4; i++) {
        Eigen::Matrix<double, -1, 1> tmp = unroll(mats[i]);
        for (int j = 0; j < tmp.rows(); j++) data[j * 4 + i] = (tmp[j]);
      }
    } else {
      // If not, I construct the dynamic projection (i.e. solve
      // systems for
      // the coeffs) and project to the superspace.

      Eigen::SparseMatrix<double> Phi =
          make_projection(knotx, knoty, xuniqueknt, yuniqueknt, xp, yp);

      for (int i = 0; i < 4; i++) {
        Eigen::Matrix<double, -1, 1> tmp =
            unroll(mats[i]).transpose() * Phi.transpose();

        for (int j = 0; j < tmp.rows(); j++) data[j * 4 + i] = (tmp[j]);
      }
    }
  }

  return;
}

/* eval() evaluates the geometry. I look up the position in the knot vector,
 * scale the input arguments, evaluate the 1D basis functions and sum over
 * them with the controll points from data. */

Eigen::Vector3d Patch::evalCore(double x, double y) const {
  const int xloc = find_loc_in_knot(x, xuniqueknt);
  const int yloc = find_loc_in_knot(y, yuniqueknt);
  const int numy = (yuniqueknt.size() - 1) * yp;
  const double scaledx = rescale(x, xuniqueknt[xloc], xuniqueknt[xloc + 1]);
  const double scaledy = rescale(y, yuniqueknt[yloc], yuniqueknt[yloc + 1]);

  double *xbasis = new double[xp];
  double *ybasis = new double[yp];

  evalBrnstnBasis(xp - 1, xbasis, scaledx);
  evalBrnstnBasis(yp - 1, ybasis, scaledy);

  double tmp[4] = {0., 0., 0., 0.};

  for (int i = 0; i < xp; i++) {
    for (int j = 0; j < yp; j++) {
      const double tpbasisval = xbasis[i] * ybasis[j];
      const int accs = 4 * (numy * (xp * xloc + i) + yp * yloc + j);
#pragma omp simd
      for (int k = 0; k < 4; k++) tmp[k] += data[accs + k] * tpbasisval;
    }
  }

  delete[] xbasis;
  delete[] ybasis;

  Eigen::Vector3d out(tmp[0], tmp[1], tmp[2]);
  // Rescaling by the NRBS weight, i.e. projection to 3D from 4D hom

  return out / tmp[3];
}

Eigen::Matrix<double, 3, 2> Patch::jacobianCore(const double x,
                                                const double y) const {
  const int xloc = find_loc_in_knot(x, xuniqueknt);
  const int yloc = find_loc_in_knot(y, yuniqueknt);
  const int numy = (yuniqueknt.size() - 1) * yp;
  const double scaledx = rescale(x, xuniqueknt[xloc], xuniqueknt[xloc + 1]);
  const double scaledy = rescale(y, yuniqueknt[yloc], yuniqueknt[yloc + 1]);

  double *xbasis = new double[xp];
  double *ybasis = new double[yp];
  double *xbasisD = new double[xp];
  double *ybasisD = new double[yp];

  evalBrnstnBasis(xp - 1, xbasis, scaledx);
  evalBrnstnBasis(yp - 1, ybasis, scaledy);
  evalBrnstnDerBasis(xp - 1, xbasisD, scaledx);
  evalBrnstnDerBasis(yp - 1, ybasisD, scaledy);

  double tmp[4] = {0., 0., 0., 0.};
  double tmpDx[4] = {0., 0., 0., 0.};
  double tmpDy[4] = {0., 0., 0., 0.};

  for (int i = 0; i < xp; i++) {
    for (int j = 0; j < yp; j++) {
      const double tpbasisval = xbasis[i] * ybasis[j];
      const double tpbasisvalDx = xbasisD[i] * ybasis[j];
      const double tpbasisvalDy = xbasis[i] * ybasisD[j];
      const int accs = 4 * (numy * (xp * xloc + i) + yp * yloc + j);

      // Here I add up the values of the basis functions in the dc
      // basis
#pragma omp simd
      for (int k = 0; k < 4; k++) {
        tmp[k] += data[accs + k] * tpbasisval;
        tmpDx[k] += data[accs + k] * tpbasisvalDx;
        tmpDy[k] += data[accs + k] * tpbasisvalDy;
      }
    }
  }

  delete[] xbasis;
  delete[] ybasis;
  delete[] xbasisD;
  delete[] ybasisD;

  Eigen::Matrix<double, 3, 2> out;

  // Eigen::Vector3d out;

  double bot = 1. / (tmp[3] * tmp[3]);

#pragma omp simd
  for (int k = 0; k < 3; k++) {
    out(k, 0) = (tmpDx[k] * tmp[3] - tmp[k] * tmpDx[3]) * bot;
    out(k, 1) = (tmpDy[k] * tmp[3] - tmp[k] * tmpDy[3]) * bot;
  }

  return out;
}

}  // namespace Spl