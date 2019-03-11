// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#include "spltest.h"

// Assumes that deBoor and Bernstein work fine.
// We initialize a random coefficient vector and act like it belongs to a
// spline. Then we do "bezier-extraction" via interpolation and check if the
// "fast" evaluation routines match the result of the deBoor recursion.

int test_localise() {
  using namespace Spl;
  constexpr int num_pts = 100;
  std::vector<double> pts;
  {
    constexpr double h = 1. / (num_pts - 1);
    for (int i = 0; i < num_pts; i++) {
      pts.push_back(i * h);
    }
  }

  for (int deg = 0; deg < 19; deg++) {
    for (int k = 0; k < 22; k++) {
      int kappa = k;
      std::vector<double> coeffs;

      {
        std::default_random_engine generator(time(0));
        std::uniform_real_distribution<double> distribution(0.0, 1.0);

        for (int i = 0; i < kappa + deg + 1; i++) {
          coeffs.push_back(distribution(generator));
        }
      }

      const auto knts = make_unif_knots(deg + 1, kappa);

      auto mask = make_interpolation_mask(deg + 1);
      const auto uniq = extract_unique_knots(knts);
      const int uniq_size = uniq.size();
      const auto pts = make_interpolation_points(uniq, mask);

      const auto vals = deBoor(coeffs, knts, pts);

      auto bezcoeffs = get_coeffs<double>(uniq_size - 1, mask, vals);

      // std::cout << "mask is ";
      // printvec(mask);
      // std::cout << "unique.size() "<< uniq_size << "\n";
      // printvec(bezcoeffs);
      // printvec(coeffs);

      const int sz = vals.size();

      std::vector<double> vals2;

      for (auto x : pts) {
        int idx = find_loc_in_knot(x, uniq);
        vals2.push_back(evalBrnstn(deg, bezcoeffs.data() + idx * (deg + 1),
                                   rescale(x, uniq[idx], uniq[idx + 1])));
      }

      for (int i = 0; i < sz; i++) {
        if (std::abs(vals[i] - vals2[i]) > 0.00000001) {
          std::cout << "oops by deg " << deg << " and kappa " << kappa
                    << std::endl;
          return 1;
        }
      }
    }
  }

  return 0;
}