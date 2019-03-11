// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef _include_knot_
#define _include_knot_

#include <vector>
namespace Spl {
/**
 *  @brief Here, routines for the creation and processing of knot vectors are
 * defined.
 *
 */

inline std::vector<double> make_bezier_knots(int deg) noexcept {
  std::vector<double> out;
  out.reserve(deg * 2);
  for (int i = 0; i < deg; i++) out.push_back(0);
  for (int i = 0; i < deg; i++) out.push_back(1);
  return out;
}

inline int get_deg_of_knts(const std::vector<double> &in) {
  constexpr double tol = .0000001;
  int i = 0;
  const int sz = in.size();
  while (++i < sz) {
    if (in[i] > tol) {
      return i;
    };
  }
  return 0;
}

inline std::vector<double> make_unif_knots(int deg, int interior,
                                           int kntrep) noexcept {
  std::vector<double> out;
  const double h = 1. / (interior + 1);
  out.reserve(deg * 2);
  for (int i = 0; i < deg; i++) out.push_back(0);
  for (int i = 1; i < interior + 1; i++)
    for (int k = 0; k < kntrep; k++) out.push_back(i * h);
  for (int i = 0; i < deg; i++) out.push_back(1);
  return out;
}

inline std::vector<double> make_unif_knots(int p, int lvl) {
  return make_unif_knots(p, lvl, 1);
}

// Chunk knots eats knotvectors and returns knotvectors in which each knot is
// unique, up to tolerance tol.
inline std::vector<double> extract_unique_knots(const std::vector<double> &in) {
  constexpr double tol = .0000001;
  std::vector<double> out;
  const int size = in.size();
  // std::cout<< "chunk " << in[0]  <<"\n";
  out.push_back(in[0]);
  for (int i = 1; i < size; i++) {
    if (in[i] > (in[i - 1] + tol)) out.push_back(in[i]);
  }
  return out;
}

// Up tp tol, fins knot element index for a given x.
inline int find_loc_in_knot(const double &x, const std::vector<double> &v) {
  constexpr double tol = .00000000001;
  int size = v.size();
  if (((x - tol) < 0.) && ((x + tol) > 0.)) return 0;
  for (int i = 0; i < size - 1; i++) {
    if (v[i] <= x && v[i + 1] > x) return i;
  }
  return size - 2;
}
}  // namespace Spl
#endif